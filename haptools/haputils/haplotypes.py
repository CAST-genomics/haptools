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

		##### Do some basic checks #####
		# hap_varids and hap_varalleles must be the same length
		assert(len(self.hap_varids)==len(self.hap_varalleles))
		# Do not repeat variant IDs
		assert(len(set(self.hap_varids))==len(self.hap_varids))
		# TODO - STR and LA checks

	def GetVariantID(self, record):
		return "%s_%s_%s_%s"%(record.CHROM, record.POS, record.REF, ":".join(record.ALT))

	def CheckHapMatch(self, vcf_reader):
		sample_list = vcf_reader.samples

		# Keep track of whether each chrom
		# of each sample matches or not
		hap_matches = {}
		for sample in sample_list:
			hap_matches[sample] = [True, True]

		found_ids = 0
		records = vcf_reader("%s:%s-%s"%(self.hap_chr, self.hap_start, self.hap_end))
		for record in records:
			# Check if the record is part of the haplotype
			if (record.ID in self.hap_varids):
				varind = self.hap_varids.index(record.ID)
			elif (GetVariantID(record) in self.varids):
				varind = self.hap_varids.index(self.GetVariantID(record))
			else: continue
			var_allele = self.hap_varalleles[varind]
			var_id = self.hap_varids[varind]
			found_ids += 1

			# Determine the index of the target allele
			record_alleles = [record.REF]+record.ALT
			assert(var_allele in record_alleles) # TODO fail more gracefully with a message
			allele_ind = record_alleles.index(var_allele)

			# Check if each sample matches
			# If not, set hap_matches to False
			for i in range(len(sample_list)):
				sample = sample_list[i]
				call = record.genotypes[i]
				assert(call[2]) # make sure it is phased
				for j in range(2):
					if call[j]!=allele_ind: hap_matches[sample][j] = False

		# Check we got through all the IDs
		assert(found_ids == len(self.hap_varids)) # TODO fail more gracefully with a message
		return hap_matches

	def Transform(self, vcf_reader):
		hap_gts = {}

		##### Case 1 - STR length #####
		pass # TODO

		##### Case 2 - variant haplotype #####
		hap_matches = self.CheckHapMatch(vcf_reader)

		##### Case 3 - LA (may be in combination with case 2) #####
		pass # TODO

		# Get haplotype counts
		for sample in vcf_reader.samples:
			hap_gts[sample] = sum(hap_matches[sample])
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