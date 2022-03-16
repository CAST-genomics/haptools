"""
Example test command:

haptools simphenotype --vcf tests/data/simple.vcf.gz --hap tests/data/simple.hap.gz \
	--out test --simu-qt  --simu-hsq 0.8 --simu-rep 2
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

	# Set up file to write summary info
	outf_sum = open("%s.par"%outprefix, "w")
	outf_sum.write("\t".join(["Haplotype", "Frequency", "Effect"])+"\n")

	# Add effect for each haplotype
	for hap in hap_reader:
		sample_to_hap = hap.Transform(vcf_reader)
		pt_sim.AddEffect(sample_to_hap, hap.hap_effect)

		# Output summary info
		hap_freq = hap.GetFrequency(vcf_reader)
		outf_sum.write("\t".join([hap.hap_id, str(hap_freq), str(hap.hap_effect)])+"\n")

	# Output simulated phenotypes to a file (haptools simphenotypes)
	pt_sim.WritePhenotypes(outprefix, simu_rep, simu_hsq, \
		simu_k, simu_qt, simu_cc)

	# Done
	outf_sum.close()