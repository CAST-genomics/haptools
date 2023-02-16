#!/usr/bin/env python

# To test: ./clumpSTR.py --summstats-snps tests/eur_gwas_pvalue_chr19.LDL.glm.linear --clump-snp-field ID --clump-field p-value --clump-chrom-field CHROM --clump-pos-field position --clump-p1 0.2 --out test.clump
import scipy.stats
import logging
import sys

from haptools.data.genotypes import GenotypesRefAlt
from haptools.data.genotypes import GenotypesTR # TODO for GenotypesTR import trtools first then if not use the code

class Variant:
    def __init__(self, varid, chrom, pos, pval, vartype):
        self.varid = varid
        self.chrom = chrom
        self.pos = pos
        self.pval = pval
        self.vartype = vartype

    def __str__(self):
        return "%s %s %s %s %s"%(self.varid, self.chrom, self.pos, self.pval, self.vartype)

class SummaryStats:
    """
    Keep track of summary statistics
    TODO: add detailed class and methods documentation
    """

    def __init__(self):
        self.summstats = []

    def Load(self, statsfile, vartype="SNP", pthresh=1.0,
            snp_field="SNP", p_field="P",
            chrom_field="CHR", pos_field="POS"):
        """
        Load summary statistics
        Ignore variants with pval < pthresh
        Not yet implemented
        """
        summstats = [] # List of Variants

        # First, parse header line to get col. numbers
        f = open(statsfile, "r")
        header_items = [item.strip() for item in f.readline().split()]
        try:
            snp_col = header_items.index(snp_field)
        except ValueError:
            print("Could not find %s in header"%snp_field)
            sys.exit(1)
        try:
            p_col = header_items.index(p_field)
        except ValueError:
            print("Could not find %s in header"%p_field)
            sys.exit(1)
        try:
            chrom_col = header_items.index(chrom_field)
        except ValueError:
            print("Could not find %s in header"%chrom_field)
            sys.exit(1)
        try:
            pos_col = header_items.index(pos_field)
        except ValueError:
            print("Could not find %s in header"%pos_field)
            sys.exit(1)

        # Now, load in stats. Skip things with pval>pthresh
        line = f.readline()
        while line.strip() != "":
            items = [item.strip() for item in line.strip().split()]
            if float(items[p_col]) > pthresh:
                line = f.readline()
                continue
            summstats.append(Variant(items[snp_col], items[chrom_col],
                int(items[pos_col]), float(items[p_col]), vartype))
            line = f.readline()
        f.close()
        self.summstats.extend(summstats)

    def GetNextIndexVariant(self, index_pval_thresh):
        """
        Get the next index variant, which is the 
        variant with the best p-value
        If no more variants below the clump-p1 threshold,
        return None
        Not yet implemented
        """
        best_var = None
        best_var_p = 1.0
        for variant in self.summstats:
            if variant.pval < best_var_p and variant.pval<index_pval_thresh:
                best_var = variant
                best_var_p = variant.pval
        return best_var

    def QueryWindow(self, indexvar, window_kb):
        """
        Find all candidate variants in the specified
        window around the index variant
        Not yet implemented
        """
        # First get stats on the indexvariant
        chrom = indexvar.chrom
        pos = indexvar.pos

        # Find candidates in the window
        candidates = []
        for variant in self.summstats:
            if variant.chrom == chrom and abs(variant.pos-pos)/1000 < window_kb:
                candidates.append(variant)

        return candidates

    def RemoveClump(self, clumpvars):
        """
        Remove the variants from a clump 
        from further consideration
        """
        keepvars = []
        for variant in self.summstats:
            if variant not in clumpvars:
                keepvars.append(variant)
        self.summstats = keepvars

def GetOverlappingSamples(snpgts, strgts):
    """
    Get overlapping set of samples
    """
    return [] # TODO

def LoadVariant(var, snpgts, strgts, samples):
    """
    Extract vector of genotypes for this variant
    """
    # if it's a SNP we should take from data matrix in GenotypesRefAlt
    # if it's a STR we should take from data matrix in GenotypesTR
    # TODO determine whether the variant is a snp or str
    return [] # TODO

def ComputeLD(candidate_gt, index_gt, snpgts, strgts, samples):
    """
    Compute the LD between two variants
    """
    # TODO - possibly check for NAs in gts and remove them
    # Compute and return Pearson r2
    return scipy.stats.pearsonr(index_gt, candidate_gt)[0]**2

def WriteClump(indexvar, clumped_vars, outf):
    """
    Write a clump to the output file
    Not yet implemented
    """
    outf.write("\t".join([indexvar.varid, indexvar.chrom, str(indexvar.pos),
        str(indexvar.pval), indexvar.vartype, 
        ",".join([str(item) for item in clumped_vars])])+"\n")


def clumpstr(summstats_snps, summstats_strs, gts_snps, gts_strs, clump_p1, clump_p2,
    clump_snp_field, clump_field, clump_chrom_field, clump_pos_field,
    clump_kb, clump_r2, out, log):
    ###### User checks ##########
    # TODO - need one of summstats_snps or summstats_strs
    # TODO - if summstats_snps, also need gts_snps
    if summstats_snps:
        assert gts_snps is not None # todo check to ensure this check works properly
    # TODO - if summstats_strs, also need gts_strs
    if summstats_strs:
        assert gts_strs is not None #TODO check to ensure this check works properly

    ###### Load summary stats ##########
    summstats = SummaryStats()
    if summstats_snps is not None:
        summstats.Load(summstats_snps, vartype="SNP", pthresh=clump_p2,
            snp_field=clump_snp_field, p_field=clump_field,
            chrom_field=clump_chrom_field, pos_field=clump_pos_field)
    if summstats_strs is not None:
        summstats.Load(summstats_strs, vartype="STR", pthresh=clump_p2,
            snp_field=clump_snp_field, p_field=clump_field,
            chrom_field=clump_chrom_field, pos_field=clump_pos_field)

    ###### Set up genotypes ##########
    snpgts = None
    strgts = None
    if gts_snps is not None:
        snpgts = GenotypesRefAlt.load(gts_snps)
    if gts_strs is not None:
        pass # TODO remove once GenotypesTR is implemented
        #strgts = GenotypesTR.load(gts_strs)
    #samples = GetOverlappingSamples(snpgts, strgts)

    ###### Setup output file ##########
    outf = open(out, "w")
    outf.write("\t".join(["ID","CHROM","POS","P","VARTYPE","CLUMPVARS"])+"\n")

    ###### Perform clumping ##########
    indexvar = summstats.GetNextIndexVariant(clump_p1)
    while indexvar is not None:
        # Load indexvar gts
        indexvar_gt = LoadVariant(indexvar, snpgts, strgts, samples)
        # Collect candidate variants within range of index variant
        candidates = summstats.QueryWindow(indexvar, clump_kb)

        # calculate LD between candidate vars and index var
        clumpvars = []
        for c in candidates:
            # load candidate variant c genotypes 
            candidate_gt = LoadVariant(c, snpgts, strgts, samplels)
            r2 = ComputeLD(candidate_gt, indexvar_gt, snpgts, strgts, samples)
            if r2 > clump_r2:
                clumpvars.append(c)
        WriteClump(indexvar, clumpvars, outf)
        summstats.RemoveClump(clumpvars+[indexvar])
        indexvar = summstats.GetNextIndexVariant(clump_p1)
    outf.close()