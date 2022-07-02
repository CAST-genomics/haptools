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
    log: Logger
        A logging instance for recording debug statements

    Examples
    --------
    >>> gens = Genotypes.load("tests/data/example.vcf.gz")
    >>> haps = Haplotypes.load("tests/data/basic.hap")
    >>> haps_gts = GenotypesRefAlt(None)
    >>> haps.transform(gens, haps_gts)
    >>> phenosim = PhenoSimulator(haps_gts)
    >>> phenotypes = phenosim.run()
    """

    def __init__(
        self,
        genotypes: Genotypes,
        output: Path = None,
        log: Logger = None,
    ):
        """
        Initialize a PhenoSimulator object

        Parameters
        ----------
        genotypes: Genotypes
            Genotypes for each haplotype
        output: Path
            Path to a '.pheno' file to which the generated phenotypes could be written
        log: Logger, optional
            A logging instance for recording debug statements
        """
        self.gens = genotypes
        self.phens = Phenotypes(fname=output)
        self.phens.names = tuple()
        self.phens.data = None
        self.phens.samples = self.gens.samples
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
        self.log.info("Computing genetic component")
        # standardize the genotypes
        gts = (gts - gts.mean(axis=0)) / gts.std(axis=0)
        # generate the genetic component
        pt = (betas * gts).sum(axis=1)
        # compute the heritability
        if heritability is None:
            self.log.info("Computing heritability")
            # if heritability is not defined, then we set it equal to the sum of the
            # effect sizes
            # assuming the genotypes are independent, this makes the variance of the
            # noise term equal to 1 - sum(betas^2)
            heritability = np.power(betas, 2).sum()
            # # account for the fact that the genotypes are not independent by adding the
            # # covariance between all of the variables
            # for a_idx, b_idx in combinations(range(len(betas)), 2):
            #     heritability += 2 * betas[a_idx] * betas[b_idx] * \
            #         np.cov(gts[:,a_idx], gts[:,b_idx])[0][1]
        self.log.info(f"Adding environmental component for h^squared: {heritability}")
        # compute the environmental effect
        noise = np.var(pt) * (np.reciprocal(heritability) - 1)
        # finally, add everything together to get the simulated phenotypes
        pt += np.random.normal(0, noise, size=pt.shape)
        if prevalence is not None:
            self.log.info(f"Converting to case/control with prevalence {prevalence}")
            # first, find the number of desired positives
            k = int(prevalence * len(pt))
            # choose the top k values and label them positive
            bool_pt = np.repeat(False, repeats=len(pt))
            bool_pt[np.argpartition(pt, k)[-k:]] = True
            pt = bool_pt
        pt = pt[:, np.newaxis]
        # now, save the archived phenotypes for later
        if self.phens.data is None:
            self.phens.data = pt
        else:
            self.phens.data = np.concatenate((self.phens.data, pt), axis=1)
        self.phens.names = self.phens.names + ("-".join(ids),)
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
