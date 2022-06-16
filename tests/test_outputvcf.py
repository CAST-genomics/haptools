import os
import pytest
import numpy as np
from cyvcf2 import VCF
from pathlib import Path
from haptools.sim_genotype import output_vcf
from haptools.admix_storage import HaplotypeSegment 
DATADIR = Path(__file__).parent.joinpath("data")

def _get_files():
    bkp_file = DATADIR.joinpath("outvcf_test.bp")
    model_file = DATADIR.joinpath("outvcf_gen.dat")
    vcf_file = DATADIR.joinpath("outvcf_test.vcf")
    sampleinfo_file = DATADIR.joinpath("outvcf_info.tab")
    out_prefix = DATADIR.joinpath("outvcf_out")
    return bkp_file, model_file, vcf_file, sampleinfo_file, out_prefix

def _get_breakpoints(bkp_file):
    # Collect breakpoints to proper format used in output_vcf function
    breakpoints = []

    with open(bkp_file, 'r') as bkp:
        sample = []
        # create breakpoints list 
        for line in bkp:
            info = line.strip().split()
            if len(info) > 1:
                # gather segments
                segment = HaplotypeSegment(info[0], int(info[1]), 
                                           int(info[2]), float(info[3]))
                sample.append(segment)
            else:
                # write out sample
                if sample:
                    breakpoints.append(sample)
                sample = []

        # end of file write out last sample
        breakpoints.append(sample)

    breakpoints = np.array(breakpoints, dtype=object)
    return breakpoints

def test_todo():
    return

def test_vcf_output():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_prefix = _get_files()
    bkps = _get_breakpoints(bkp_file)

    # generate output vcf file
    output_vcf(bkps, model_file, vcf_file, sampleinfo_file, str(out_prefix))

    # Expected output for each variant (note these are phased so order matters)
    # CHROM	POS  FORMAT	 Sample1      Sample2
    # 1	10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1	59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2	10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    # read in vcf file
    vcf = VCF(str(out_prefix)+'.vcf')
    for var in vcf:
        if var.CHROM == "1" and var.POS == 10114:
            assert(var.genotypes[0] == [0,0,True])
            assert(var.format('POP')[0] == "YRI,YRI")
            assert(var.genotypes[1] == [1,1,True])
            assert(var.format('POP')[1] == "CEU,CEU")

        elif var.CHROM == "1" and var.POS == 59423090:
            assert(var.genotypes[0] == [0,1,True])
            assert(var.format('POP')[0] == "CEU,YRI")
            assert(var.genotypes[1] == [1,0,True])
            assert(var.format('POP')[1] == "YRI,CEU")

        elif var.CHROM == "2" and var.POS == 10122:
            assert(var.genotypes[0] == [1,0,True])
            assert(var.format('POP')[0] == "YRI,CEU")
            assert(var.genotypes[1] == [0,1,True])
            assert(var.format('POP')[1] == "CEU,YRI")

        else:
            assert(False)

    # Remove output file from output_vcf located at out_prefix + '.vcf'
    os.remove(str(out_prefix)+'.vcf')
    return

def test_mismatched_chroms():
    # Test to ensure that the code errors when a chromosomes is found in one of the vcf or bp file but not the other.
    return

