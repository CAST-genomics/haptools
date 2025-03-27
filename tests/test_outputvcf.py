import os
from pathlib import Path

import pytest
import numpy as np
from cyvcf2 import VCF

from haptools.logging import getLogger
from haptools.data import GenotypesPLINK
from haptools.admix_storage import HaplotypeSegment
from haptools.sim_genotype import (
    _prepare_coords,
    output_vcf,
    validate_params,
    simulate_gt,
    write_breakpoints,
)

DATADIR = Path(__file__).parent.joinpath("data")


def _get_files(plink_input=False, plink_output=False):
    log = getLogger(name="test")
    bkp_file = DATADIR / "outvcf_test.bp"
    model_file = DATADIR / "outvcf_gen.dat"
    vcf_file = DATADIR / ("outvcf_test" + (".pgen" if plink_input else ".vcf.gz"))
    sampleinfo_file = DATADIR / "outvcf_info.tab"
    out_file = DATADIR / ("outvcf_out" + (".pgen" if plink_output else ".vcf.gz"))
    return bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log


def _get_random_files():
    log = getLogger(name="test")
    model_file = DATADIR / "outvcf_gen_random.dat"
    vcf_file = DATADIR / "outvcf_test_random.vcf.gz"
    sampleinfo_file = DATADIR / "outvcf_info_random.tab"
    out_file = DATADIR / "outvcf_out_random.vcf.gz"
    out_prefix = DATADIR / "outvcf_out_random"
    coords_dir = DATADIR / "map"
    return model_file, vcf_file, sampleinfo_file, coords_dir, out_prefix, out_file, log


def _get_breakpoints(bkp_file, model_file):
    # Collect breakpoints to proper format used in output_vcf function
    breakpoints = []
    # create pop_dict
    mfile = open(model_file, "r")
    num_samples, *pops = mfile.readline().strip().split()
    num_samples = int(num_samples)
    pop_dict = {}
    for pop_ind, pop in enumerate(pops):
        pop_dict[pop] = pop_ind

    with open(bkp_file, "r") as bkp:
        sample = []
        # create breakpoints list
        for line in bkp:
            info = line.strip().split()
            if len(info) > 1:
                # gather segments
                segment = HaplotypeSegment(
                    pop_dict[info[0]], int(info[1]), int(info[2]), float(info[3])
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


def _get_expected_output():
    gts = GenotypesPLINK(None)
    # CHROM POS  FORMAT  Sample_1      Sample_2
    # 1 10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1 59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2 10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    gts.data = np.zeros(12).reshape((3, 2, 2)).astype(np.uint8)
    gts.data[[1, 2], [1, 0], 0] = 1
    gts.data[[1, 2], [0, 1], 1] = 1
    gts.data[0, 1, :] = 1
    gts.data = gts.data.transpose((1, 0, 2))
    gts.variants = np.array(
        [
            ("1:10114:T:C", "1", 10114, ("T", "C")),
            ("1:59423090:A:G", "1", 59423090, ("A", "G")),
            ("2:10122:A:G", "2", 10122, ("A", "G")),
        ],
        dtype=gts.variants.dtype,
    )
    gts.samples = ("Sample_1", "Sample_2")
    return gts


def test_end_bkp_coords():
    coords_dir = DATADIR / "map"
    chroms = ["22"]
    region = False
    coords, np_coords, max_coords, end_coords = _prepare_coords(
        coords_dir, chroms, region
    )
    assert coords[0][-1].get_bp_pos() == np.iinfo(np.int32).max
    assert end_coords[0].get_bp_pos() == np.iinfo(np.int32).max


def test_variants_greater_than_last_coord():
    log = getLogger(name="test")
    bkp_file = DATADIR / "var_greater.bp"
    vcf_file = DATADIR / "var_greater.vcf.gz"
    model_file = DATADIR / "outvcf_gen.dat"
    sampleinfo_file = DATADIR / "outvcf_info.tab"
    out_file = DATADIR / "outvcf_out.vcf.gz"
    chroms = ["1"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output vcf file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        True,
        False,
        str(out_file),
        log,
    )

    # read in vcf file
    vcf = VCF(str(out_file))
    for var in vcf:
        if var.CHROM == "1" and var.POS == 10114:
            assert var.genotypes[0] == [0, 0, True]
            assert var.format("POP")[0] == "YRI,YRI"
            assert var.genotypes[1] == [1, 1, True]
            assert var.format("POP")[1] == "CEU,CEU"
        elif var.CHROM == "1" and var.POS == 249403765:
            assert var.genotypes[0] == [0, 1, True]
            assert var.format("POP")[0] == "CEU,YRI"
            assert var.genotypes[1] == [1, 0, True]
            assert var.format("POP")[1] == "YRI,CEU"
        else:
            assert False

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_alt_chrom_name():
    # Test when the ref VCF has chr{X|\d+} form
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    bkp_file = DATADIR / "outvcf_test_chr.bp"
    vcf_file = DATADIR / "outvcf_test_chr.vcf"
    chroms = ["1", "2", "X"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output vcf file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        True,
        False,
        str(out_file),
        log,
    )

    # read in vcf file
    vcf = VCF(str(out_file))
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

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_no_replace():
    # Test too few samples to generate a VCF when sampling without replacement
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    bkp_file = DATADIR / "outvcf_test_no_replace.bp"
    vcf_file = DATADIR / "outvcf_test_no_replace.vcf"
    chroms = ["1", "2", "X"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # test when we have enough CEU/YRI samples to simulate
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        True,
        True,
        str(out_file),
        log,
    )

    # read in vcf file
    vcf = VCF(str(out_file))
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

    # generate output vcf file should fail due to lack of CEU samples
    # this info file has only one CEU individual
    with pytest.raises(Exception, match="No available sample"):
        one_ceu_sampleinfo_file = DATADIR / "outvcf_info-one_CEU.tab"
        bkp_file = DATADIR / "outvcf_test_no_replace2.bp"
        bkps = _get_breakpoints(bkp_file, model_file)
        output_vcf(
            bkps,
            chroms,
            model_file,
            str(vcf_file),
            one_ceu_sampleinfo_file,
            None,
            True,
            True,
            True,
            str(out_file),
            log,
        )

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_vcf_output():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    chroms = ["1", "2"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output vcf file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        True,
        False,
        str(out_file),
        log,
    )

    # Expected output for each variant (note these are phased so order matters)
    # CHROM	POS  FORMAT	 Sample1      Sample2
    # 1	10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1	59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2	10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    # read in vcf file
    vcf = VCF(str(out_file))
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

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_vcf_randomness():
    chroms = ["1"]
    region = {
        "chr": "1",
        "start": int("1"),
        "end": int("10115"),
    }
    (
        model_file,
        vcf_file,
        sampleinfo_file,
        coords_dir,
        out_prefix,
        out_file,
        log,
    ) = _get_random_files()

    # create random breakpoints file with around 1000 output samples to ensure we get all genotypes output
    num_samples, pop_dict, breakpoints = simulate_gt(
        model_file,
        coords_dir,
        chroms,
        region,
        10000,
        log,
    )
    bkps = write_breakpoints(num_samples, pop_dict, breakpoints, str(out_prefix), log)
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        False,
        False,
        str(out_file),
        log,
    )

    # If random haplotypes arent chosen per person we would only have an output gt of 0|1
    #     for same population simulated individuals
    vcf = VCF(str(out_file))
    for var in vcf:
        all_gts = set()
        for gt, sample_pops in zip(var.genotypes, var.format("POP")):
            sample_pops = sample_pops.split(",")
            if sample_pops[0] == sample_pops[1]:
                all_gts.add(np.sum(gt[:2]))
    assert sorted(list(all_gts)) == [0, 1, 2]

    # Clean up by removing the output file from output_vcf
    out_prefix.with_suffix(".bp").unlink()
    out_file.unlink()


def test_someflags_vcf():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    chroms = ["1", "2"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output vcf file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        True,
        False,
        False,
        str(out_file),
        log,
    )

    # Expected output for each variant (note these are phased so order matters)
    # CHROM POS  FORMAT  Sample1      Sample2
    # 1 10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1 59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2 10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    # read in vcf file
    vcf = VCF(str(out_file))
    for var in vcf:
        if var.CHROM == "1" and var.POS == 10114:
            assert var.genotypes[0] == [0, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT
            assert var.genotypes[1] == [1, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT

        elif var.CHROM == "1" and var.POS == 59423090:
            assert var.genotypes[0] == [0, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT
            assert var.genotypes[1] == [1, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT

        elif var.CHROM == "2" and var.POS == 10122:
            assert var.genotypes[0] == [1, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT
            assert var.genotypes[1] == [0, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" in var.FORMAT

        else:
            assert False

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_noflags_vcf():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    chroms = ["1", "2"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output vcf file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        False,
        False,
        False,
        str(out_file),
        log,
    )

    # Expected output for each variant (note these are phased so order matters)
    # CHROM POS  FORMAT  Sample1      Sample2
    # 1 10114    GT:POP  0|0:YRI,YRI  1|1:CEU,CEU
    # 1 59423090 GT:POP  0|1:CEU,YRI  1|0:YRI,CEU
    # 2 10122    GT:POP  1|0:YRI,CEU  0|1:CEU,YRI
    # read in vcf file
    vcf = VCF(str(out_file))
    for var in vcf:
        if var.CHROM == "1" and var.POS == 10114:
            assert var.genotypes[0] == [0, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT
            assert var.genotypes[1] == [1, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT

        elif var.CHROM == "1" and var.POS == 59423090:
            assert var.genotypes[0] == [0, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT
            assert var.genotypes[1] == [1, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT

        elif var.CHROM == "2" and var.POS == 10122:
            assert var.genotypes[0] == [1, 0, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT
            assert var.genotypes[1] == [0, 1, True]
            assert "SAMPLE" not in var.FORMAT
            assert "POP" not in var.FORMAT

        else:
            assert False

    # Clean up by removing the output file from output_vcf
    out_file.unlink()


def test_pgen_output():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files(
        plink_output=True
    )
    chroms = ["1", "2"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output pgen file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        False,
        False,
        False,
        str(out_file),
        log,
    )

    expected = _get_expected_output()
    gts = GenotypesPLINK(out_file)
    gts.read()
    gts.check_phase()

    np.testing.assert_allclose(gts.data, expected.data)
    assert gts.samples == expected.samples
    assert np.array_equal(gts.variants, expected.variants)

    out_file.unlink()
    out_file.with_suffix(".pvar").unlink()
    out_file.with_suffix(".psam").unlink()


def test_pgen_input():
    # read in all files and breakpoints
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files(
        plink_input=True, plink_output=True
    )
    chroms = ["1", "2"]
    bkps = _get_breakpoints(bkp_file, model_file)

    # generate output pgen file
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        None,
        False,
        False,
        False,
        str(out_file),
        log,
    )

    expected = _get_expected_output()
    gts = GenotypesPLINK(out_file)
    gts.read()
    gts.check_phase()

    np.testing.assert_allclose(gts.data, expected.data)
    assert gts.samples == expected.samples
    assert np.array_equal(gts.variants, expected.variants)

    out_file.unlink()
    out_file.with_suffix(".pvar").unlink()
    out_file.with_suffix(".psam").unlink()


def test_region_bkp():
    modelfile = DATADIR / "outvcf_gen.dat"
    popsize = 100000
    region = {"chr": "22", "start": 16000, "end": 18000}
    coords_dir = DATADIR / "map"
    chroms = ["22"]
    seed = 100
    log = getLogger(name="test")
    num_samples, pop_dict, all_samples = simulate_gt(
        modelfile, coords_dir, chroms, region, popsize, log, seed
    )

    # Make sure lowest bkp listed is 16111 and greatest is 18674
    for sample in all_samples:
        for coord in sample:
            assert (16111 <= coord.get_end_coord() <= 18674) or (
                coord.get_end_coord() == np.iinfo(np.int32).max
            )


def test_region_vcf():
    region = {"chr": "2", "start": 1, "end": 10122}
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files()
    bkps = _get_breakpoints(bkp_file, model_file)
    chroms = ["2"]
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        region,
        True,
        True,
        False,
        str(out_file),
        log,
    )

    vcf = VCF(str(out_file))
    for var in vcf:
        assert var.POS == 10122 and var.CHROM == "2"
        assert var.genotypes[0] == [1, 0, True]
        assert var.format("POP")[0] == "YRI,CEU"
        assert var.genotypes[1] == [0, 1, True]
        assert var.format("POP")[1] == "CEU,YRI"

    out_file.unlink()


def test_region_pgen():
    region = {"chr": "2", "start": 1, "end": 10122}
    bkp_file, model_file, vcf_file, sampleinfo_file, out_file, log = _get_files(
        plink_input=True
    )
    bkps = _get_breakpoints(bkp_file, model_file)
    chroms = ["2"]
    output_vcf(
        bkps,
        chroms,
        model_file,
        str(vcf_file),
        sampleinfo_file,
        region,
        False,
        False,
        False,
        str(out_file),
        log,
    )

    vcf = VCF(str(out_file))
    for var in vcf:
        assert var.POS == 10122 and var.CHROM == "2"
        assert var.genotypes[0] == [1, 0, True]
        assert var.genotypes[1] == [0, 1, True]

    out_file.unlink()


# model_file exception validation
def test_model_files():
    mapdir = DATADIR / "map"
    chroms = ["50"]

    popsize = 1000
    vcf_file = str(DATADIR / "outvcf_test.vcf")
    sampleinfo_file = DATADIR / "outvcf_info.tab"

    faulty_model = DATADIR / "dat_files/faulty_model_sample_number_to_int.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Can't convert samples number to an integer."

    faulty_model = DATADIR / "dat_files/faulty_model_num_pops.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (
        str(e.value)
    ) == "Invalid number of populations given: 1. We require at least 2."

    faulty_model = DATADIR / "dat_files/faulty_model_less_than_1.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Number of samples is less than 1."

    # validate exception number of pops = number in pop_fracs
    faulty_model = DATADIR / "dat_files/faulty_model_pop_fracs.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Can't convert generation to integer."

    faulty_model = DATADIR / "dat_files/faulty_model_frac.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Can't convert population fractions to type float."

    faulty_model = DATADIR / "dat_files/faulty_model_pop_header.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (
        (str(e.value))
        == "Total fractions given to populations do not match number of populations in"
        " the header."
    )

    faulty_model = DATADIR / "dat_files/faulty_model_cur_gen.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (
        (str(e.value)) == "Current generation 1 - previous generation 4 = -3"
        " is less than 1. Please ensure the generations given in the first column"
        " are correct."
    )

    faulty_model = DATADIR / "dat_files/faulty_model_sum_frac.dat"
    with pytest.raises(Exception) as e:
        validate_params(
            faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Population fractions for generation 1 do not sum to 1."

    # Validate mapdir exceptions
    model = DATADIR / "dat_files/correct_model.dat"
    faulty_mapdir = DATADIR / "maps"
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Map directory given is not a valid path."

    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == f"Chromosome {chroms[0]} in the list given is not valid."

    chroms = ["1"]
    faulty_mapdir = DATADIR / "test_map"
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (
        (str(e.value))
        == "No valid coordinate files found. Must contain chr{1-22,X} in the file name"
        " and end in .map"
    )

    faulty_mapdir = DATADIR / "test_map_2"
    with pytest.raises(Exception) as e:
        validate_params(
            model, faulty_mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
        )
    assert (
        (str(e.value))
        == "No valid coordinate files found. Must contain chr{1-22,X} in the file name"
        " and end in .map"
    )

    # validate popsize exceptions
    faulty_popsize = "NA"
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, faulty_popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Popsize is not an Integer."

    faulty_popsize = 0
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, faulty_popsize, vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Popsize must be greater than 0."

    # validate vcf sample collection exception
    faulty_vcf_file = DATADIR / "faulty_vcf.vcf"
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, faulty_vcf_file, sampleinfo_file, False
        )
    assert (str(e.value)) == "Unable to collect vcf samples."

    # validate sample_info file exception
    faulty_sampleinfo_file = DATADIR / "faulty_info.tab"
    with pytest.raises(Exception) as e:
        validate_params(
            model, mapdir, chroms, popsize, vcf_file, faulty_sampleinfo_file, False
        )
    msg = (
        "Sample HG00022 from population CEU in sampleinfo file is not present in the"
        " vcf file."
    )
    assert (str(e.value)) == msg

    faulty_model = DATADIR / "dat_files/faulty_model_sample_info.dat"
    mfile = open(faulty_model, "r")
    num_samples, *pops = mfile.readline().strip().split()
    with pytest.raises(Exception) as e:
        for model_pop in pops:
            validate_params(
                faulty_model, mapdir, chroms, popsize, vcf_file, sampleinfo_file, False
            )
        msg = (
            f"Population {model_pop} in model file is not present in the sample info"
            " file."
        )
        assert str(e.value) == msg
