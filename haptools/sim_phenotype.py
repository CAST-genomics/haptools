from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass, field
from logging import getLogger, Logger, DEBUG

import numpy as np
import numpy.typing as npt

from .data import Haplotype as HaplotypeBase
from .data import (
    Extra,
    Genotypes,
    Phenotypes,
    Haplotypes,
    GenotypesPLINK,
    GenotypesRefAlt,
)


@dataclass
class Haplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for simphenotype

    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    beta: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(Extra("beta", ".2f", "Effect size in linear model"),),
    )


class PhenoSimulator:
    """
    Simulate phenotypes from genotypes

    Attributes
    ----------
    gens: Genotypes
        Genotypes to simulate
    phens: Phenotypes
        Simulated phenotypes; filled by :py:meth:`~.PhenoSimular.run`
    rng: np.random.Generator, optional
        A numpy random number generator
    log: Logger
        A logging instance for recording debug statements

    Examples
    --------
    >>> gens = Genotypes.load("tests/data/example.vcf.gz")
    >>> haps = Haplotypes.load("tests/data/basic.hap")
    >>> haps_gts = haps.transform(gens)
    >>> phenosim = PhenoSimulator(haps_gts)
    >>> phenosim.run(haps.data.values())
    >>> phenotypes = phenosim.phens
    """

    def __init__(
        self,
        genotypes: Genotypes,
        output: Path = Path("/dev/stdout"),
        seed: int = None,
        log: Logger = None,
    ):
        """
        Initialize a PhenoSimulator object

        Parameters
        ----------
        genotypes: Genotypes
            Genotypes for each haplotype
        output: Path, optional
            Path to a '.pheno' file to which the generated phenotypes could be written

            Defaults to stdout if not provided
        seed: int, optional
            A seed to initialize the random number generator

            This is useful if you want the generated phenotypes to be the same across
            multiple PhenoSimulator instances. If not provided, it will be random.
        log: Logger, optional
            A logging instance for recording debug statements
        """
        self.gens = genotypes
        self.phens = Phenotypes(fname=output)
        self.phens.data = None
        self.phens.samples = self.gens.samples
        self.rng = np.random.default_rng(seed)
        self.log = log or getLogger(self.__class__.__name__)

    def run(
        self,
        effects: list[Haplotype],
        heritability: float = None,
        prevalence: float = None,
    ) -> npt.NDArray:
        """
        Simulate phenotypes for an entry in the Genotypes object

        The generated phenotypes will also be added to
        :py:attr:`~.PhenoSimulator.output`

        Parameters
        ----------
        effects: list[Haplotype]
            A list of Haplotypes to use in an additive fashion within the simulations
        heritability: float, optional
            The simulated heritability of the trait

            If not provided, this will default to the sum of the squared effect sizes
        prevalence: float, optional
            How common should the disease be within the population?

            If this value is specified, case/control phenotypes will be generated
            instead of quantitative traits.

        Returns
        -------
        npt.NDArray
            The simulated phenotypes, as a np array of shape num_samples x 1
        """
        # extract the relevant haplotype info from the Haplotype objects
        ids = [hap.id for hap in effects]
        betas = np.array([hap.beta for hap in effects])
        self.log.debug(f"Extracting haplotype genotypes for haps: {ids}")
        self.log.debug(f"Beta values are {betas}")
        # extract the haplotype "genotypes" and compute the phenotypes
        gts = self.gens.subset(variants=ids).data[:, :, :2].sum(axis=2)
        self.log.info(f"Computing genetic component w/ {gts.shape[0]} causal effects")
        # standardize the genotypes
        std = gts.std(axis=0)
        gts = (gts - gts.mean(axis=0)) / std
        # for genotypes where the stdev is 0, just set all values to 0 instead of nan
        zero_elements = std == 0
        gts[:, zero_elements] = np.zeros((gts.shape[0], np.sum(zero_elements)))
        # generate the genetic component
        pt = (betas * gts).sum(axis=1)
        # compute the heritability
        if heritability is None:
            self.log.debug("Computing heritability as the sum of the squared betas")
            heritability = np.power(betas, 2).sum()
            if heritability > 1:
                heritability = 1
            # compute the environmental effect
            noise = 1 - heritability
        else:
            # compute the environmental effect
            noise = np.var(pt) * (np.reciprocal(heritability) - 1)
        self.log.info(f"Adding environmental component {noise} for h^2 {heritability}")
        # finally, add everything together to get the simulated phenotypes
        pt_noise = self.rng.normal(0, np.sqrt(noise), size=pt.shape)
        if self.log.getEffectiveLevel() == DEBUG:
            # if we're in debug mode, compute the pearson correlation and report it
            # but don't do this otherwise to keep things fast
            corr = np.corrcoef(pt, pt + pt_noise)[1, 0]
            self.log.debug(f"Estimated heritability is {corr}")
        pt += pt_noise
        # now, handle case/control
        name_suffix = ""
        if prevalence is not None:
            self.log.info(f"Converting to case/control with prevalence {prevalence}")
            # first, find the number of desired positives
            k = int(prevalence * len(pt))
            # choose the top k values and label them positive
            bool_pt = np.zeros(len(pt), dtype=np.bool_)
            if k == len(pt):
                bool_pt = ~bool_pt
            else:
                max_indices = np.argpartition(-pt, k)[:k]
                bool_pt[max_indices] = True
            pt = bool_pt
            name_suffix = "-cc"
        # now, save the archived phenotypes for later
        self.phens.append(name="-".join(ids) + name_suffix, data=pt.astype(np.float64))
        return pt

    def write(self):
        """
        Write the generated phenotypes to the file specified in
        :py:meth:`~.PhenoSimular.__init__`
        """
        self.phens.write()


def simulate_pt(
    genotypes: Path,
    haplotypes: Path,
    num_replications: int = 1,
    heritability: float = None,
    prevalence: float = None,
    region: str = None,
    samples: list[str] = None,
    haplotype_ids: set[str] = None,
    chunk_size: int = None,
    output: Path = Path("-"),
    log: Logger = None,
):
    """
    Haplotype-aware phenotype simulation. Create a set of simulated phenotypes from a
    set of haplotypes.

    GENOTYPES must be formatted as a VCF or PGEN file and HAPLOTYPES must be formatted
    according to the .hap format spec

    Note: GENOTYPES must be the output from the the transform subcommand.

    \f
    Examples
    --------
    >>> haptools simphenotype tests/data/example.vcf.gz tests/data/example.hap.gz > simu_phens.tsv

    Parameters
    ----------
    genotypes : Path
        The path to the transformed genotypes in VCF or PGEN format
    haplotypes : Path
        The path to the haplotypes in a .hap file
    replications : int, optional
        The number of rounds of simulation to perform
    heritability : int, optional
        The heritability of the simulated trait; must be a float between 0 and 1

        If not provided, it will be computed from the sum of the squared effect sizes
    prevalence : int, optional
        The prevalence of the disease if the trait should be simulated as case/control;
        must be a float between 0 and 1

        If not provided, a quantitative trait will be simulated, instead
    region : str, optional
        The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'

        For this to work, the VCF and .hap file must be indexed and the seqname must
        match!

        Defaults to loading all haplotypes
    sample : tuple[str], optional
        A subset of the samples from which to extract genotypes

        Defaults to loading genotypes from all samples
    samples_file : Path, optional
        A single column txt file containing a list of the samples (one per line) to
        subset from the genotypes file
    haplotype_ids: set[str], optional
        A list of haplotype IDs to obtain from the .hap file. All others are ignored.

        If not provided, all haplotypes will be used.
    chunk_size: int, optional
        The max number of variants to fetch from the PGEN file at any given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory. This argument is ignored if the genotypes are
        not in PGEN format.
    output : Path, optional
        The location to which to write the simulated phenotypes
    log : Logger, optional
        The logging module for this task
    """
    if log is None:
        log = logging.getLogger("haptools simphenotype")
        logging.basicConfig(
            format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
            level="ERROR",
        )

    log.info("Loading haplotypes")
    hp = Haplotypes(haplotypes, haplotype=Haplotype, log=log)
    hp.read(region=region, haplotypes=haplotype_ids)

    if haplotype_ids is None:
        haplotype_ids = set(hp.data.keys())

    if genotypes.suffix == ".pgen":
        log.info("Loading genotypes from PGEN file")
        gt = GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
    else:
        log.info("Loading genotypes from VCF/BCF file")
        gt = Genotypes(fname=genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=haplotype_ids)
    log.info("QC-ing genotypes")
    gt.check_missing()
    gt.check_biallelic()

    # check that all of the genotypes were loaded successfully and warn otherwise
    if len(haplotype_ids) < len(gt.variants):
        diff = list(haplotype_ids.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} haplotypes could not be found in the genotypes file. Check "
            "that the hap IDs in your .hap file correspond with those in the genotypes"
            f" file. Here are the first few missing variants: {diff[:first_few]}"
        )

    # Initialize phenotype simulator (haptools simphenotype)
    log.info("Simulating phenotypes")
    pt_sim = PhenoSimulator(gt, output=output, log=log)
    for i in range(num_replications):
        pt_sim.run(hp.data.values(), heritability, prevalence)
    log.info("Writing phenotypes")
    pt_sim.write()
