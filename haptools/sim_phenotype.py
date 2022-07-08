from __future__ import annotations
from pathlib import Path
from itertools import combinations
from logging import getLogger, Logger
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from .data import Haplotype as HaplotypeBase
from .data import GenotypesRefAlt, Phenotypes, Haplotypes, Extra


@dataclass
class Haplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for simphenotype

    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    ancestry: str
    beta: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(
            Extra("ancestry", "s", "Local ancestry"),
            Extra("beta", ".2f", "Effect size in linear model"),
        ),
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
        heritability: float = 1,
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

            If not provided, this will be estimated from the variability of the
            genotypes
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
        self.log.info(f"Extracting haplotype genotypes for haps: {ids}")
        # extract the haplotype "genotypes" and compute the phenotypes
        gts = self.gens.subset(variants=ids).data.sum(axis=2)
        self.log.info(f"Computing genetic component w/ {gts.shape[0]} causal effects")
        # standardize the genotypes
        gts = (gts - gts.mean(axis=0)) / gts.std(axis=0)
        # generate the genetic component
        pt = (betas * gts).sum(axis=1)
        # compute the heritability
        self.log.info(f"Adding environmental component for h^2: {heritability}")
        # compute the environmental effect
        noise = np.var(pt) * (np.reciprocal(heritability) - 1)
        # finally, add everything together to get the simulated phenotypes
        pt += self.rng.normal(0, noise, size=pt.shape)
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
    heritability: float = 1,
    prevalence: float = None,
    region: str = None,
    samples: list[str] = None,
    output: Path = Path("-"),
    log: Logger = None,
):
    """
    Haplotype-aware phenotype simulation. Create a set of simulated phenotypes from a
    set of haplotypes.

    GENOTYPES must be formatted as a VCF and HAPLOTYPES must be formatted according
    to the .hap format spec

    \f
    Examples
    --------
    >>> haptools simphenotype tests/data/example.vcf.gz tests/data/example.hap.gz > simu_phens.tsv

    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF format
    haplotypes : Path
        The path to the haplotypes in a .hap file
    replications : int, optional
        The number of rounds of simulation to perform
    heritability : int, optional
        The heritability of the simulated trait; must be a float between 0 and 1
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
    output : Path, optional
        The location to which to write the simulated phenotypes
    log : Logger, optional
        The logging module for this task
    """
    if log is None:
        log = logging.getLogger("run")
        logging.basicConfig(
            format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
            level="ERROR",
        )

    log.info("Loading haplotypes")
    hp = Haplotypes(haplotypes, haplotype=Haplotype, log=log)
    hp.read(region=region)

    log.info("Extracting variants from haplotypes")
    variants = {var.id for hap in hp.data.values() for var in hap.variants}

    log.info("Loading genotypes")
    gt = GenotypesRefAlt(genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    log.info("QC-ing genotypes")
    gt.check_missing()
    gt.check_biallelic()
    gt.check_phase()

    log.info("Transforming genotypes via haplotypes")
    hp_gt = GenotypesRefAlt(fname=None, log=log)
    hp.transform(gt, hp_gt)

    # Initialize phenotype simulator (haptools simphenotype)
    log.info("Simulating phenotypes")
    pt_sim = PhenoSimulator(hp_gt, output=output, log=log)
    for i in range(num_replications):
        pt_sim.run(hp.data.values(), heritability, prevalence)
    log.info("Writing phenotypes")
    pt_sim.write()
