#!/usr/bin/env python

# To test: ./clumpSTR.py --summstats-snps tests/eur_gwas_pvalue_chr19.LDL.glm.linear --clump-snp-field ID --clump-field p-value --clump-chrom-field CHROM --clump-pos-field position --clump-p1 0.2 --out test.clump
import scipy.stats
import numpy as np
import logging
import math
import sys

from haptools.data.genotypes import GenotypesRefAlt
from haptools.data.genotypes import GenotypesTR # TODO for GenotypesTR import trtools first then if not use the code

# TODO update all print statements with log variable

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
        header = f.readline()
        if header.startswith('#'):
            header = header[1:]
        header_items = [item.strip() for item in header.split()]
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

def LoadVariant(var, snpgts, strgts, log):
    """
    Extract vector of genotypes for this variant
    """
    # if it's a SNP we should take from data matrix in GenotypesRefAlt
    # if it's a STR we should take from data matrix in GenotypesTR
    # Grab variant from snps or strs depending on variant type
    if var.vartype.lower() == 'snp':
        var_ind = (snpgts.variants['pos']==int(var.pos)) & \
                  (snpgts.variants['chrom']==var.chrom)
        variant_gts = np.sum(snpgts.data[:, var_ind,:], axis=2).flatten()
    else:
        var_ind = (strgts.variants['pos']==int(var.pos)) & \
                  (strgts.variants['chrom']==var.chrom)
        variant_gts = np.sum(strgts.data[:, var_ind,:], axis=2).flatten()
    
    return variant_gts

def _CalcChiSQ(f00, f01, f10, f11, gt_counts, n):
    """
    Calculate Chi-squared test stat for given freqs.
    """
    chisq_exp = np.zeros((3,3))
    root_exp = np.zeros((3,3))

    # calculate expected values for a given root
    root_exp[0,0] = n * f00**2
    root_exp[0,1] = 2 * n * f00 * f01
    root_exp[0,2] = n * f01**2
    root_exp[1,0] = 2 * n * f00 * f10
    root_exp[1,1] = 2 * n * f01 * f10 + 2 * n * f00 * f11 
    root_exp[1,2] = 2 * n * f01 * f11
    root_exp[2,0] = n * f10**2
    root_exp[2,1] = 2 * n * f10 * f11
    root_exp[2,2] = n * f11**2
    for i in range(3):
        for j in range(3):
            if root_exp[i,j] > 0.0:
                chisq_exp = ((gt_counts[i,j] - root_exp[i,j])**2/root_exp[i,j])

    return np.sum(chisq_exp)

def _CalcLDStats(f00, p, q, gt_counts, n):
    """
    Given frequency of gt 0|0 (f11) and major and minor allele freqs p and q calculate stats.
    """
    f01 = p - f00
    f10 = q - f00
    f11 = 1 - (f00 + f01 + f10)
    D = (f00 * f11) - (f01 * f10)
    if D >= 0.0:
        Dmax = min(p*(1.0-q), q*(1.0-p))
    else:
        Dmax = min(p*q,(1-p)*(1-q))
    Dprime = D/Dmax
    r_squared = (D**2)/(p*(1-p)*q*(1-q))

    return round(Dprime,6), \
           round(r_squared,6), \
           _CalcChiSQ(f00, f01, f10, f11, gt_counts, n)

def _CalcBestRoot(real_roots, minhap, maxhap, p, q, gt_counts, n):
    """
    Given a list of real roots (max possible 3) calculate the best root.
    The best root is the one with the lowest chisq test statistic value.
    """
    # determine the best root by grabbing the one with the lowest chisq
    best_Dprime = 0
    best_rsquared = 0
    best_chisq = np.inf
    for root_freq00 in real_roots:
        # calculate LD stats given root freq is within bounds
        if root_freq00 >= minhap - 0.00001 and root_freq00 <= maxhap + 0.00001:
            Dprime, r_squared, chisq = _CalcLDStats(root_freq00, p, q, gt_counts, n)
            if chisq < best_chisq:
                best_Dprime = Dprime
                best_rsquared = r_squared
                best_chisq = chisq
        
    return best_Dprime, best_rsquared

def ComputeMlsLD(candidate_gt, index_gt, log):
    """
    Compute maximum likelihood solution of haplotype frequencies to calculate r squared value.
    # TODO add https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-428 to docs
    # TODO add https://github.com/t0mrg/cubex to docs

    NOTE currently this approach only works for biallelic variants since having more variants
    causes the equation we're solving for to be a cubic but instead to the degree of n where 
    n is the total number of alelles which also invalidates STRs.

    Parameters
    ----------
    candidate_gt: np.array
        array of size (genotypes,) where genotypes is the number of samples
    index_gt: np.array
        array of size (genotypes,) where genotypes is the number of samples

    Returns
    -------
    r_squared: float
        R squared value inferred from ML solution.
    """
    # load in 3x3 array where axes are genotypes (0,1,2) for each variant
    # y-axis = candidate gt, x-axis = index_gt
    gt_counts = np.zeros((3,3))
    for gt1 in range(3):
        # subset candidate gts to genotype gt
        subset_gt1 = candidate_gt == gt1
        for gt2 in range(3):
            subset_index_gt = index_gt[subset_gt1]
            gt_counts[gt1,gt2] = np.sum(subset_index_gt == gt2)

    n = np.sum(gt_counts)
    p = (2.0*np.sum(gt_counts[0,:]) + np.sum(gt_counts[1,:]))/(2.0 * n)
    q = (2.0*np.sum(gt_counts[:,0]) + np.sum(gt_counts[:,1]))/(2.0 * n)

    num_alt = (2.0*gt_counts[0,0] + gt_counts[0,1] + gt_counts[1,0])
    a = 4.0*n
    b = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*num_alt - gt_counts[1,1]
    c = -num_alt*(1.0 - 2.0*p - 2.0*q) - gt_counts[1,1]*(1.0 - p - q) + 2.0*n*p*q
    d = -num_alt*p*q

    minhap = num_alt / (2.0 * float(n))
    maxhap = (num_alt + gt_counts[1,1]) / (2.0 * float(n))
    
    xN = -b/(3.0*a)
    d2 = (math.pow(b,2)-3.0*a*c)/(9*math.pow(a,2))
    yN = a * math.pow(xN,3) + b * math.pow(xN,2) + c * xN + d
    yN2 = math.pow(yN,2)
    h2 = 4 * math.pow(a,2) * math.pow(d2,3)

    # store all real roots to cubic to iterate over and determine which is best
    real_roots = []

    # three possible scenarios of solutions
    if yN2 > h2:
        # calculate real root alpha
        number1 = 0.0
        number2 = 0.0
        if (1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))) < 0:
            number1 = -math.pow(-(1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: 
            number1 = math.pow((1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        
        if (1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))) < 0:
            number2 = -math.pow(-(1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: 
            number2 = math.pow((1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)

        # singular real root
        alpha = xN + number1 + number2

        # store real root alpha
        real_roots = [alpha]

    elif yN2 == h2:
        # Calculate three real roots alpha beta and gamma
        delta = math.pow((yN/2.0*a),(1.0/3.0))
        alpha = xN + delta
        beta = xN + delta
        gamma = xN - 2.0*delta

        # store all real roots
        real_roots = [alpha, beta, gamma]

    elif yN2 < h2:
        # calculate 3 real roots alpha beta and gamma
        h = math.pow(h2, 0.5)
        theta = ((math.acos(-yN/h))/3.0)
        delta = math.pow(d2,0.5)
        alpha = xN + 2.0 * delta * math.cos(theta)
        beta = xN + 2.0 * delta * math.cos(2.0 * math.pi/3.0 + theta)
        gamma = xN + 2.0 * delta * math.cos(4.0 * math.pi/3.0 + theta)

        # store all real roots
        real_roots = [alpha, beta, gamma]

    else: 
        raise Exception(f"Can't calculate r squared from given values {yN2} and {h2}")

    # Solve for best roots
    best_Dprime, best_rsquared = _CalcBestRoot(real_roots, minhap, maxhap, p, q, gt_counts, n)
    return best_Dprime, best_rsquared

def ComputeLD(candidate_gt, index_gt, LD_type, log):
    """
    Compute the LD between two variants
    """
    # TODO - Check for NAs in gts and remove them
    # Compute and Maximum likelihood solution or Pearson r2
    if LD_type == 'MLS':
        return ComputeMlsLD(candidate_gt, index_gt, log)
    elif LD_type == 'Pearson':
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
    clump_kb, clump_r2, LD_type, out, log):
    ###### User checks ##########
    # TODO NEED TO ADD THE LOGGER TO EACH FUNCTION TO TRACK DEBUG AND INFO MESSAGES
    # TODO - need one of summstats_snps or summstats_strs
    # TODO - if summstats_snps, also need gts_snps
    if summstats_snps:
        assert gts_snps is not None # TODO check to ensure this check works properly
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
        log.debug("Loading SNP Genotypes.")
        snpgts = GenotypesRefAlt.load(gts_snps) # TODO need to check if input is PGEN (use GenotypesPLINK instead)
    if gts_strs is not None:
        log.debug("Loading STR Genotypes.")
        # TODO Implement STR loader 
        #strgts = GenotypesTR.load(gts_strs)
    if gts_snps and gts_strs:
        log.debug("Calculating set of overlapping samples between STRs and SNPs.")
        samples = GetOverlappingSamples(snpgts, strgts)

    # NOTE snpgts has data, variants, and samples where data is alleles (samples x variants x alleles)
    #      variants has id, pos, chrom, ref, alt for snps
    #      samples is list of samples corresponding to x-axis of samples
    # TODO uncomment subsample snps and strs to set of overlapping samples
    #snpgts.data = snpgts.data[samples,:,:]
    #snpgts.samples = snpgts.samples[samples]
    #strgts.data = snpgts.data[samples,:,:]
    #strgts.samples = snpgts.samples[samples]

    ###### Setup output file ##########
    outf = open(out, "w")
    outf.write("\t".join(["ID","CHROM","POS","P","VARTYPE","CLUMPVARS"])+"\n")

    ###### Perform clumping ##########
    indexvar = summstats.GetNextIndexVariant(clump_p1)
    while indexvar is not None:
        # Load indexvar gts
        indexvar_gt = LoadVariant(indexvar, snpgts, strgts, log)
        # Collect candidate variants within range of index variant
        candidates = summstats.QueryWindow(indexvar, clump_kb)

        log.debug(f"Current index variant: {indexvar}")

        # calculate LD between candidate vars and index var
        clumpvars = []
        for c in candidates:
            # load candidate variant c genotypes 
            candidate_gt = LoadVariant(c, snpgts, strgts, log)
            Dprime, r2 = ComputeLD(candidate_gt, indexvar_gt, LD_type, log)
            log.debug(
                    f"D\' and r2 between {indexvar} with {c}\n" +
                    f"D\' = {Dprime}, r^2 = {r2}"
            )
            if r2 > clump_r2:
                clumpvars.append(c)
        WriteClump(indexvar, clumpvars, outf)
        summstats.RemoveClump(clumpvars+[indexvar])
        indexvar = summstats.GetNextIndexVariant(clump_p1)
    outf.close()