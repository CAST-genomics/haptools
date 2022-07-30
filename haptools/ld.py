from __future__ import annotations
import logging
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from haptools import data
from .data import Haplotype as HaplotypeBase


@dataclass
class Haplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for the ld command

    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    ld: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(data.Extra("ld", ".3f", "Linkage-disequilibrium"),),
    )


def calc_ld(
    target: str,
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: list[str] = None,
    haplotype_ids: set[str] = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    output: Path = Path("/dev/stdout"),
    log: Logger = None,
):
    """
    Creates a VCF composed of haplotypes

    Parameters
    ----------
    target : str
        The ID of the haplotype or variant with which we will calculate LD
    genotypes : Path
        The path to the genotypes
    haplotypes : Path
        The path to the haplotypes in a .hap file
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
        and :py:meth:`~.data.Haplotypes.read`
    samples : list[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    haplotype_ids: set[str], optional
        A set of haplotype IDs to obtain from the .hap file. All others are ignored.

        If not provided, all haplotypes will be used.
    chunk_size: int, optional
        The max number of variants to fetch from the PGEN file at any given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory. This argument is ignored if the genotypes are
        not in PGEN format.
    discard_missing : bool, optional
        Discard any samples that are missing any of the required genotypes

        The default is simply to complain about it
    output : Path, optional
        The location to which to write output
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = logging.getLogger("haptools ld")
        logging.basicConfig(
            format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
            level="ERROR",
        )

    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes, log=log)
    if haplotype_ids is not None:
        haplotype_ids.add(target)
    hp.read(region=region, haplotypes=haplotype_ids)

    log.info("Extracting variants from haplotypes")
    variants = {var.id for hap in hp.data.values() for var in hap.variants}

    # check to see whether the target was a haplotype
    if target in hp.data:
        target = hp.data.pop(target)
        if len(hp.data) == 0:
            log.error(
                "There must be at least one more haplotype in the .hap file "
                "than the TARGET haplotype specified."
            )

    # check that all of the haplotypes were loaded successfully and warn otherwise
    if haplotype_ids is not None and len(haplotype_ids) > len(hp.data):
        diff = list(haplotype_ids.difference(hp.data.keys()))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} haplotypes could not be found in the .hap file. Check "
            "that the IDs in your .hap file correspond with those you provided. "
            f"Here are the first few missing haplotypes: {diff[:first_few]}"
        )

    # if the target was not a haplotype, then it must be a variant
    # so we must load it from the genotype file
    if not isinstance(target, data.Haplotype):
        variants.add(target)

    if genotypes.suffix == ".pgen":
        log.info("Loading genotypes from PGEN file")
        gt = data.GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
    else:
        log.info("Loading genotypes from VCF/BCF file")
        gt = data.GenotypesRefAlt(fname=genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing(discard_also=discard_missing)
    gt.check_biallelic()
    gt.check_phase()

    # check that all of the variants were loaded successfully and warn otherwise
    if len(variants) < len(gt.variants):
        # check to see whether the target got loaded if it wasn't a haplotype
        if not (isinstance(target, data.Haplotype) or (target in gt.variants["id"])):
            raise ValueError(
                "Could not find the provided target ID among either the haplotypes "
                "in the .hap file or the variants in the genotype file. Check that "
                f"'{target}' appears in either the .hap file or the genotype file."
            )
        # report the missing variants
        diff = list(variants.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} variants could not be found in the genotypes file. Check "
            "that the IDs in your .hap file correspond with those in the genotypes "
            f"file. Here are the first few missing variants: {diff[:first_few]}"
        )

    log.info("Transforming genotypes via haplotypes")
    hp_gt = data.GenotypesRefAlt(fname=None, log=log)
    hp.transform(gt, hp_gt)

    log.info("Obtaining target genotypes")
    if isinstance(target, data.Haplotype):
        target_gts = target.transform(gt).sum(axis=1)
    else:
        target_gts = gt.subset(variants=(target,)).data[:, 0, :2].sum(axis=1)

    log.info("Computing LD between haplotypes and the target")
    # construct a new Haplotypes object that also stores the LD values
    hp_out = data.Haplotypes(fname=output, haplotype=Haplotype, log=log)
    hp_out.data = {}
    for hap_id in hp.data:
        # break the BaseHaplotype instance up into its properties
        hapd = hp.data[hap_id].__dict__
        hapd_variants = hapd.pop("variants")
        # obtain genotypes for the current haplotype
        hap_gts = hp_gt.subset(variants=(hap_id,)).data[:, 0, :2].sum(axis=1)
        # compute the LD between the current haplotype and the target
        hapd["ld"] = np.corrcoef(target_gts, hap_gts)[1, 0]
        # create the Haplotype instance and add the variants in
        hp_out.data[hap_id] = Haplotype(**hapd)
        hp_out.data[hap_id].variants = hapd_variants

    log.info("Outputting .hap file with LD values")
    hp_out.write()
