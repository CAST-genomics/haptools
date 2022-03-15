"""
Classes for reading/writing and other functions
on haplotypes

TODO:
* transform
  * deal with LA and STRs

* writephenotypes
* documentation and type hinting
* tests
"""
import gzip
import tabix

class Haplotype:
	def __init__(self, hap_id, hap_chr, hap_start, hap_end, \
					hap_varids, hap_varalleles, \
					hap_LA, hap_effect):
		self.hap_id = hap_id
		self.hap_chr = hap_chr
		self.hap_start = hap_start
		self.hap_end = hap_end
		self.hap_varids = hap_varids
		self.hap_varalleles = hap_varalleles
		self.hap_LA = hap_LA
		self.hap_effect = hap_effect

	def Transform(self, vcf_reader):
		hap_gts = {}
		for sample in vcf_reader.samples:
			hap_gts[sample] = 0 # TODO sample -> 0/1/2 for the haplotype
		return hap_gts

	def __str__(self):
		return self.hap_id

class HapReader:
	def __init__(self, hapfile, region=None):
		if region is None:
			self.haps = gzip.open(hapfile)
		else:
			self.haps = tabix.open(hapfile).querys(region)

	def GetHaplotype(self, hapline):
		hap_id = hapline[0]
		hap_chr = hapline[1]
		hap_start = int(hapline[2])
		hap_end = int(hapline[3])
		hap_varids = hapline[4].split(",")
		hap_varalleles = hapline[5].split(",")
		hap_LA = hapline[6].split(",")
		hap_effect = float(hapline[7])
		return Haplotype(hap_id, hap_chr, hap_start, hap_end, \
							hap_varids, hap_varalleles, \
							hap_LA, hap_effect)

	def __iter__(self):
		return self

	def __next__(self):
		# Get the next haplotype
		hapline = next(self.haps)

		# If reading from gzip file, convert to list
		# so it is in the same format as from tabix
		if type(hapline) == bytes:
			hapline = hapline.decode('utf-8')
			while hapline.startswith('#'):
				hapline = next(self.haps).decode('utf-8')
			hapline = hapline.strip().split()

		# Conver the line to a haplotype object
		return self.GetHaplotype(hapline)

class PhenoSimulator:
	def __init__(self, sample_list):
		self.pts = {}
		for sample in sample_list:
			self.pts[sample] = 0

	def AddEffect(self, sample_to_hap, effect):
		for sample in self.pts.keys():
			assert(sample in sample_to_hap.keys())
			self.pts[sample] += effect*sample_to_hap[sample]

	def WritePhenotypes(self, outprefix, simu_rep, simu_hsq, \
						simu_k, simu_qt, simu_cc):
		pass # TODO write output files