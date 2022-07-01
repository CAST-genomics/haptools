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
        log: Logger = None,
    ):
        """
        Initialize a PhenoSimulator object

        Parameters
        ----------
        genotypes: Genotypes
            Genotypes for each haplotype
        log: Logger, optional
            A logging instance for recording debug statements
        """
        self.gens = genotypes
        self.log = log or getLogger(self.__class__.__name__)

    def run(
        self,
        varID: str,
        beta: float,
        heritability: float = 1,
        prevalence: float = None,
    ) -> npt.NDArray:
        """
        Simulate phenotypes for an entry in the Genotypes object

        Parameters
        ----------
        varID: str
            The ID of the variant in the Genotype object from which to simulate
            phenotypes
        beta: float
            The simulated effect size of the variant; must be between -1 and 1
        heritability: float, optional
            The simulated heritability of the trait
        prevalence: float, optional
            How common should the disease be within the population?

            If this value is specified, case/control phenotypes will be generated
            instead of quantitative traits.
        """
        gts = self.gens.subset(variants=[varID]).data[:, 0].sum(axis=1)
        # error =
        pt = effect * (gts - gts.mean(axis=0)) / gts.std(axis=0)

        # Simulate and write
        for i in range(simu_rep):
            resid_component = np.random.normal(0, np.var(pt) * (1 / simu_hsq - 1))
            if simu_qt:
                outinfo.append(pt + resid_component)
            else:
                raise NotImplementedError("Case control not implemented")


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
    pt_sim = PhenoSimulator(hp_gt)
    phens = pt_sim.run()

    # Set up file to write summary info
    outf_sum = open("%s.par" % outprefix, "w")
    outf_sum.write("\t".join(["Haplotype", "Frequency", "Effect"]) + "\n")

    # Add effect for each haplotype
    for idx, hap in enumerate(haps.data.values()):
        sample_to_hap = hap.transform(gens)
        pt_sim.add_effect(sample_to_hap, hap.beta)

        # Output summary info
        hap_freq = hap.GetFrequency(gens)
        outf_sum.write("\t".join([hap.id, str(hap_freq), str(hap.beta)]) + "\n")

    # Output simulated phenotypes to a file (haptools simphenotypes)
    pt_sim.run(outprefix, simu_rep, simu_hsq, simu_k, simu_qt, simu_cc)

    # Done
    outf_sum.close()
