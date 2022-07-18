import os
import pytest
import numpy as np
from cyvcf2 import VCF
from pathlib import Path
from haptools.sim_genotype import output_vcf, validate_params
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

class TestValidateParams():

    #model files
    #read it in as open(model_file, 'r')
 
    

    #model_file exception validation
    def test_model_files(self):
        faulty_model = DATADIR.joinpath("dat_files/faulty_model.dat")
        mapdir = DATADIR.joinpath("plink.chr1.GRCh38.map")
        chroms = 1
        popsize = 0
        vcf_file = DATADIR.joinpath("outvcf_test.vcf")
        sampleinfo_file = DATADIR.joinpath("outvcf_info.tab")

        with pytest.raises(Exception) as e:
            validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
            assert (str(e.value)) == "Can't convert samples number to an integer."
            assert (str(e.value)) == "Invalid number of populations given: {num_pops}. We require at least 2."

     
        faulty_model = DATADIR.joinpath("dat_files/faul_mod.dat")
        with pytest.raises(Exception):
            validate_params(faulty_model, mapdir, chroms, .1, vcf_file, sampleinfo_file)
            assert Exception == "Number of samples is less than 1."
        

   
"""""
    

    
     with pytest.raises(Exception):
        validate_params(faulty_model, mapdir, chroms, .1, vcf_file, sampleinfo_file)
        assert Exception == "Number of samples is less than 1."

        

  
          #when code fails it throws an exception -- make sure proper exception is thrown
    model_file.read()
        #sum of each row must equal 1 
    if model_file > 1 or model_file < 1:
        assert(False)

        #admixed first like must = 0
    assert ["Admixed"][1] == 0

        #each row must be greater than 
        #code that handels this gives exception
        #need to refer to the specific column i am counting
    count_row = 0
    for count in model_file:
        if count - count_row != 0:
            assert(false)
        count_row + 1


    #maps
    #we care about 1,3,4 column
    #just make sure they exist

    #output
    #make sure it has some sort of prefix - path + string at end of path
    #this is a path
    ##simple - 2 column (for now)
    #population code




        try:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    except Exception as e:
        assert e == Exception("Can't convert samples number to an integer.")


    
"""

    
