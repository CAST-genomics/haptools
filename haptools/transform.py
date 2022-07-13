from __future__ import annotations
import logging
from pathlib import Path

from haptools import data


def transform_haps(
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: list[str] = None,
    haplotype_ids: set[str] = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    output: Path = Path("-"),
    log: Logger = None,
):
    """
    Creates a VCF composed of haplotypes

    Parameters
    ----------
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
        Discard any samples that are missing any of the required samples

        The default is simply to complain about it
    output : Path, optional
        The location to which to write output
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = logging.getLogger("run")
        logging.basicConfig(
            format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
            level="ERROR",
        )

    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes, log=log)
    hp.read(region=region, haplotypes=haplotype_ids)

    log.info("Extracting variants from haplotypes")
    variants = {var.id for hap in hp.data.values() for var in hap.variants}

    if genotypes.suffix == ".pgen":
        log.info("Loading genotypes from PGEN file")
        gt = data.GenotypesPLINK(genotypes, log=log, chunk_size=chunk_size)
    else:
        log.info("Loading genotypes from VCF/BCF file")
        gt = data.GenotypesaRefAlt(genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing(discard_also=discard_missing)
    gt.check_biallelic()
    gt.check_phase()

    # check that all of the variants were loaded successfully and warn otherwise
    if len(variants) < len(gt.variants):
        diff = list(variants.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} variants could not be found in the genotypes file. Check "
            "that the IDs in your .hap file correspond with those in the genotypes "
            f"file. Here are the first few missing variants: {diff[:first_few]}"
        )

    log.info("Transforming genotypes via haplotypes")
    hp_gt = data.GenotypesRefAlt(fname=output, log=log)
    hp.transform(gt, hp_gt)

    log.info("Writing haplotypes to VCF")
    hp_gt.write()

    return hp_gt
