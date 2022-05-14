from __future__ import annotations
from logging import getLogger, Logger

from haplotype import HaptoolsHaplotype as Haplotype
from .data import GenotypesRefAlt, Phenotypes, Haplotypes


class PhenoSimulator:
    """
    Simulate phenotypes from genotypes

    Attributes
    ----------
    gens: GenotypesRefAlt
        A set of genotypes to use when simulating phenotypes
    haps: Haplotypes
        A set of haplotypes to use as causal variables in the simulation
    haps_gts: GenotypesRefAlt
        Genotypes for each haplotype
    log: Logger
        A logging instance for recording debug statements
    """

    def __init__(
        self, genotypes: GenotypesRefAlt, haplotypes: Haplotypes, log: Logger = None
    ):
        self.gens = genotypes
        self.haps = haplotypes
        self.log = log or getLogger(self.__class__.__name__)
        self.haps_gts = self.haps.transform(self.gens)

    def run(
        self,
        simu_rep: int,
        simu_hsq: float,
        simu_k: float,
        simu_qt: bool = False,
        simu_cc: bool = False,
    ) -> Phenotypes:
        effect = np.array([hap.beta for hap in haps.data.values()])
        gts = self.haps_gts.data.sum(axis=2)
        pt = np.sum(effect * (gts - gts.mean(axis=0)) / gts.std(axis=0))
        phens = np.empty((len(self.genotypes.samples)), dtype=np.float64)

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
    phenotypes: Path,
    simu_rep: int,
    simu_hsq: float,
    simu_k: float,
    simu_qt: bool = False,
    simu_cc: bool = False,
):
    # Initialize readers
    gens = GenotypesRefAlt.load(genotypes)
    haps = Haplotypes(haplotypes, haplotype=Haplotype)
    haps.read()

    # Initialize phenotype simulator (haptools simphenotypes)
    pt_sim = PhenoSimulator(gens, haps)
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
