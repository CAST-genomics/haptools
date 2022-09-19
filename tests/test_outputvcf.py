import os
import pytest
import numpy as np
from cyvcf2 import VCF
from pathlib import Path
from haptools.sim_genotype import output_vcf, validate_params, simulate_gt
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

    with open(bkp_file, "r") as bkp:
        sample = []
        # create breakpoints list
        for line in bkp:
            info = line.strip().split()
            if len(info) > 1:
                # gather segments
                segment = HaplotypeSegment(
                    info[0], int(info[1]), int(info[2]), float(info[3])
                )
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

def test_alt_chrom_name():
    # Test when the ref VCF has chr{X|\d+} form
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_prefix = _get_files()
    bkp_file = DATADIR.joinpath("outvcf_test_chr.bp")
    vcf_file = DATADIR.joinpath("outvcf_test_chr.vcf")
    chroms = ['1', '2', 'X']
    bkps = _get_breakpoints(bkp_file)

    # generate output vcf file
    output_vcf(bkps, chroms, model_file, vcf_file, sampleinfo_file, None, str(out_prefix))

    # read in vcf file
    vcf = VCF(str(out_prefix) + ".vcf")
    for var in vcf:
        if var.CHROM == "chr1" and var.POS == 10114:
            assert var.genotypes[0] == [0, 0, True]
            assert var.format("POP")[0] == "YRI,YRI"
            assert var.genotypes[1] == [1, 1, True]
            assert var.format("POP")[1] == "CEU,CEU"

        elif var.CHROM == "chr1" and var.POS == 59423090:
            assert var.genotypes[0] == [0, 1, True]
            assert var.format("POP")[0] == "CEU,YRI"
            assert var.genotypes[1] == [1, 0, True]
            assert var.format("POP")[1] == "YRI,CEU"

        elif var.CHROM == "chr2" and var.POS == 10122:
            assert var.genotypes[0] == [1, 0, True]
            assert var.format("POP")[0] == "YRI,CEU"
            assert var.genotypes[1] == [0, 1, True]
            assert var.format("POP")[1] == "CEU,YRI"

        elif var.CHROM == "chrX" and var.POS == 10122:
            assert var.genotypes[0] == [1, 1, True]
            assert var.format("POP")[0] == "YRI,YRI"
            assert var.genotypes[1] == [1, 1, True]
            assert var.format("POP")[1] == "YRI,YRI"

        else:
            assert False

    # Remove output file from output_vcf located at out_prefix + '.vcf'
    os.remove(str(out_prefix) + ".vcf")
    return

def test_vcf_output():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_prefix = _get_files()
    chroms = ['1', '2']
    bkps = _get_breakpoints(bkp_file)

    # generate output vcf file
    output_vcf(bkps, chroms, model_file, vcf_file, sampleinfo_file, None, str(out_prefix))

    # Expected output for each variant (note these are phased so order matters)
    # CHROM	POS  FORMAT	 Sample1      Sample2
    # 1	10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1	59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2	10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    # read in vcf file
    vcf = VCF(str(out_prefix) + ".vcf")
    for var in vcf:
        if var.CHROM == "1" and var.POS == 10114:
            assert var.genotypes[0] == [0, 0, True]
            assert var.format("POP")[0] == "YRI,YRI"
            assert var.genotypes[1] == [1, 1, True]
            assert var.format("POP")[1] == "CEU,CEU"

        elif var.CHROM == "1" and var.POS == 59423090:
            assert var.genotypes[0] == [0, 1, True]
            assert var.format("POP")[0] == "CEU,YRI"
            assert var.genotypes[1] == [1, 0, True]
            assert var.format("POP")[1] == "YRI,CEU"

        elif var.CHROM == "2" and var.POS == 10122:
            assert var.genotypes[0] == [1, 0, True]
            assert var.format("POP")[0] == "YRI,CEU"
            assert var.genotypes[1] == [0, 1, True]
            assert var.format("POP")[1] == "CEU,YRI"

        else:
            assert False

    # Remove output file from output_vcf located at out_prefix + '.vcf'
    os.remove(str(out_prefix) + ".vcf")
    return

def test_region_bkp():
    modelfile = DATADIR.joinpath("outvcf_gen.dat")
    popsize = 100000
    region = {'chr':'22','start':16000, 'end':18000}
    coords_dir = DATADIR.joinpath("map")
    chroms = ["22"]
    seed = 100
    num_samples, all_samples = simulate_gt(modelfile, coords_dir, chroms, region, popsize, seed)
    
    # Make sure lowest bkp listed is 16111 and greatest is 18674
    for sample in all_samples:
        for coord in sample:
            assert 16111 <= coord.get_end_coord() <= 18674
    return

def test_region_vcf():
    region = {'chr':'2', 'start':1, 'end':10122}
    bkp_file, model_file, vcf_file, sampleinfo_file, out_prefix = _get_files()
    bkps = _get_breakpoints(bkp_file)
    chroms = ['2']
    output_vcf(bkps, chroms, model_file, vcf_file, sampleinfo_file, region, str(out_prefix))

    vcf = VCF(str(out_prefix) + ".vcf")
    for var in vcf:
        assert var.POS == 10122 and var.CHROM == '2'
        assert var.genotypes[0] == [1, 0, True]
        assert var.format("POP")[0] == "YRI,CEU"
        assert var.genotypes[1] == [0, 1, True]
        assert var.format("POP")[1] == "CEU,YRI"

    os.remove(str(out_prefix) + ".vcf")
    return

# model_file exception validation
def test_model_files():
    mapdir = DATADIR.joinpath("map")
    chroms = ["50"]

    popsize = 1000
    vcf_file = DATADIR.joinpath("outvcf_test.vcf")
    sampleinfo_file = DATADIR.joinpath("outvcf_info.tab")

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_sample_number_to_int.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Can't convert samples number to an integer."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_num_pops.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (
        str(e.value)
    ) == "Invalid number of populations given: 1. We require at least 2."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_less_than_1.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Number of samples is less than 1."

    # validate exception number of pops = number in pop_fracs
    faulty_model = DATADIR.joinpath("dat_files/faulty_model_pop_fracs.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Can't convert generation to integer."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_frac.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Can't convert population fractions to type float."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_pop_header.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (
        (str(e.value))
        == "Total fractions given to populations do not match number of populations in"
        " the header."
    )

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_cur_gen.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (
        (str(e.value))
        == "Current generation 1 - previous generation 4 = -3"
        " is less than 1. Please ensure the generations given in the first column"
        " are correct."
    )

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_sum_frac.dat")
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (
        str(e.value)
    ) == "Population fractions for generation 1 do not sum to 1."

    # Validate mapdir exceptions
    model = DATADIR.joinpath("dat_files/correct_model.dat")
    faulty_mapdir = DATADIR.joinpath("maps")
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Map directory given is not a valid path."

    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == f"Chromosome {chroms[0]} in the list given is not valid."

    chroms = ["1"]
    faulty_mapdir = DATADIR.joinpath("test_map")
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Could not parse map directory files."

    faulty_mapdir = DATADIR.joinpath("test_map_2")
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
        )
    assert (
        (str(e.value))
        == "No valid coordinate files found. Must contain chr{1-22,X} in the file"
        " name."
    )

    # validate popsize exceptions
    faulty_popsize = "NA"
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, faulty_popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Popsize is not an Integer."

    faulty_popsize = 0
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, faulty_popsize, vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Popsize must be greater than 0."

    # validate vcf sample collection exception
    faulty_vcf_file = DATADIR.joinpath("faulty_vcf.vcf")
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, faulty_vcf_file, sampleinfo_file, None
        )
    assert (str(e.value)) == "Unable to collect vcf samples."

    # validate sample_info file exception
    faulty_sampleinfo_file = DATADIR.joinpath("faulty_info.tab")
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, vcf_file, faulty_sampleinfo_file, None
        )
    assert (str(e.value)) == "Sample HG00022 in sampleinfo file is not present in the vcf file."

    faulty_model = DATADIR.joinpath("dat_files/faulty_model_sample_info.dat")
    mfile = open(faulty_model, "r")
    num_samples, *pops = mfile.readline().strip().split()
    with pytest.raises(Exception) as e:
        for model_pop in pops:
            validate_params(
                faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, None
            )
        assert (
            (str(e.value))
            == f"Population {model_pop} in model file is not present in the sample info"
            " file."
        )
