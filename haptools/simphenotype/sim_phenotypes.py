"""
Example test command:
haptools simphenotype --vcf tests/data/simple.vcf.gz \
  --hap tests/data/simple.hap.gz \
  --out test
"""

from cyvcf2 import VCF
from ..haputils.haplotypes import *

def simulate_pt(vcffile, hapfile, simu_rep, \
        simu_hsq, simu_k, simu_qt, simu_cc, outprefix):

	# Initialize readers
	vcf_reader = VCF(vcffile, gts012=True)
	hap_reader = HapReader(hapfile)

	# Initialize phenotype simulator (haptools simphenotypes)
	pt_sim = PhenoSimulator(vcf_reader.samples)

	# Add effect for each haplotype
	for hap in hap_reader:
		sample_to_hap = hap.Transform(vcf_reader)
		pt_sim.AddEffect(sample_to_hap, hap.hap_effect)

	# Output simulated phenotypes to a file (haptools simphenotypes)
	pt_sim.WritePhenotypes(outprefix, simu_rep, simu_hsq, \
		simu_k, simu_qt, simu_cc)
