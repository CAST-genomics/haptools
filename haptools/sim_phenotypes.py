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
        # extract the ID and effect size information from the Haplotype objects
        ids = [hap.id for hap in effects]
        betas = np.array([hap.beta for hap in effects])
        # extract the haplotype "genotypes" and compute the phenotypes
        gts = self.gens.subset(variants=ids).data.sum(axis=2)
        # standardize the genotypes
        gts = (gts - gts.mean(axis=0)) / gts.std(axis=0)
        # generate the genetic component
        pt = (betas * gts).sum(axis=1)
        # compute the heritability
        if heritability is None:
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
        # compute the environmental effect
        noise = np.var(pt) * (np.reciprocal(heritability) - 1)
        # finally, add everything together to get the simulated phenotypes
        pt += np.random.normal(0, noise, size=pt.shape)
        # TODO: implement case/control phenotypes
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
    simu_rep: int,
    simu_hsq: float,
    simu_k: float,
    region: str = None,
    samples: list[str] = None,
    output: Path = Path("-"),
    log: Logger = None,
):
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
    pt_sim.run(hp.data.values())
    log.info("Writing phenotypes")
    pt_sim.write()
