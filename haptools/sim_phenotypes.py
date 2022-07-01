from __future__ import annotations
from logging import getLogger, Logger

import numpy.typing as npt

from haplotype import HaptoolsHaplotype as Haplotype
from .data import GenotypesRefAlt, Phenotypes, Haplotypes


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
        varID: str,
        beta: float,
        heritability: float = None,
        prevalence: float = None,
    ) -> npt.NDArray:
        """
        Simulate phenotypes for an entry in the Genotypes object

        The generated phenotypes will also be added to :py:attr:`~.PhenoSimulator.fname`

        Parameters
        ----------
        varID: str
            The ID of the variant in the Genotype object from which to simulate
            phenotypes
        beta: float
            The simulated effect size of the variant; must be between -1 and 1
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
        """
        gts = self.gens.subset(variants=[varID]).data[:, 0].sum(axis=1)
        pt = effect * (gts - gts.mean(axis=0)) / gts.std(axis=0)
        # add an error term
        if heritability:
            # as defined in GTCA
            pt += np.random.normal(0, np.var(pt) * ((1/heritability) - 1) )
        else:
            pt += np.random.normal(0, 1 - np.var(gts)**2)
        if self.phens.data is None:
            self.phens.data = pt
        else:
            self.phens.data = np.concatenate((self.phens.data, pt), axis=1)
        self.phens.names = self.phens.names + (varID,)
        # TODO: implement case/control phenotypes
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
    hp = data.Haplotypes(haplotypes, log=log)
    hp.read(region=region)

    log.info("Extracting variants from haplotypes")
    variants = {var.id for hap in hp.data.values() for var in hap.variants}

    log.info("Loading genotypes")
    gt = data.GenotypesRefAlt(genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing()
    gt.check_biallelic()
    gt.check_phase()

    log.info("Transforming genotypes via haplotypes")
    hp_gt = data.GenotypesRefAlt(fname=None, log=log)
    hp.transform(gt, hp_gt)

    # Initialize phenotype simulator (haptools simphenotype)
    log.info("Simulating phenotypes")
    pt_sim = PhenoSimulator(hp_gt, output=output, log=log)
    for hap in hp.data.values():
        pt_sim.run(hap.id, hap.beta)
    pt_sim.write()
