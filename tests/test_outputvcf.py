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



#model_file exception validation
def test_model_files():
    mapdir = DATADIR.joinpath("map/")
    print(mapdir)
    chroms = 1
    popsize = 1000
    vcf_file = DATADIR.joinpath("outvcf_test.vcf")
    sampleinfo_file = DATADIR.joinpath("outvcf_info.tab")

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_sample_number_to_int.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Can't convert samples number to an integer."


    faulty_model = DATADIR.joinpath("dat_files/faulty_model_num_pops.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Invalid number of populations given: {num_pops}. We require at least 2."


    faulty_model = DATADIR.joinpath("dat_files/faulty_model_less_than_1.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Number of samples is less than 1."

#validate number of pops = number in pop_fracs
    faulty_model = DATADIR.joinpath("dat_files/faulty_model_pop_fracs.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Can't convert generation to integer."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_frac.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Can't convert population fractions to type float."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_pop_header.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Total fractions given to populations do not match number of populations in the header."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_cur_gen.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Current generation {cur_gen} - previous generation {prev_gen} = {sim_gens} is less than 1. ""Please ensure the generations given in the first column are correct."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_sum_frac.dat")
    with pytest.raises(Exception) as e:
        validate_params(faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Population fractions for generation {cur_gen} do not sum to 1."

    # Validate mapdir ensuring it contains proper files.
    model = DATADIR.joinpath("dat_files/haiti.dat")
    faulty_mapdir = DATADIR.joinpath("maps")
    with pytest.raises(Exception) as e:
        validate_params(model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Map directory given is not a valid path."

    faulty_mapdir = DATADIR.joinpath("test_map/")
    with pytest.raises(Exception) as e:
        validate_params(model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Could not parse map directory files."

    faulty_mapdir = DATADIR.joinpath("test_map_2/")
    with pytest.raises(Exception) as e:
        validate_params(model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "No valid coordinate files found. Must contain chr\{1-22,X\} in the file name."

    
    popsize = "NA"
    with pytest.raises(Exception) as e:
        validate_params(model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == "Popsize is not an Integer."



    """""
    chroms = 50
    with pytest.raises(Exception) as e:
        validate_params(model, mapdir, chroms, popsize, vcf_file, sampleinfo_file)
    assert (str(e.value)) == f"Chromosome {chroms} in the list given is not valid."
    """

