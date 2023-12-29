from pathlib import Path

import pytest
import numpy as np
from cyvcf2 import VCF
from click.testing import CliRunner

from haptools.__main__ import main
from haptools.logging import getLogger
from haptools.data.genotypes import GenotypesVCF
from haptools.data.genotypes import GenotypesTR
from haptools.clump import (
    GetOverlappingSamples,
    _SortSamples,
    LoadVariant,
    _FilterGts,
    ComputeLD,
    clumpstr,
    Variant,
)

DATADIR = Path(__file__).parent.joinpath("data")
log = getLogger(name="test")


class TestClump:
    def _ld_expected(self):
        return [0.5625, 0.6071, 0.5977, 0.5398]

    def test_loading_snps(self):
        gts_snps = DATADIR / "outvcf_test.vcf.gz"
        snpgts = GenotypesVCF.load(str(gts_snps))
        snpvars = [
            Variant("test1", "1", "10114", "0.05", "snp"),
            Variant("test2", "1", "59423090", "0.05", "snp"),
            Variant("test3", "2", "10122", "0.05", "snp"),
        ]
        strgts = None

        answers = [
            np.array([2, 2, 0, 0, 0]),
            np.array([0, 0, 2, 2, 2]),
            np.array([0, 0, 2, 2, 2]),
        ]

        for var, answer in zip(snpvars, answers):
            vargts = LoadVariant(var, snpgts, log)
            vargts = np.sum(vargts, axis=1).flatten()
            assert len(vargts) == snpgts.data.shape[0]
            assert np.array_equal(vargts, answer)

    def test_sample_sorting(self):
        test1 = ["Sample_02", "Sample_00", "Sample_01"]
        test2 = ["Sample_3", "Sample_2", "Sample_1"]
        test3 = ["Sample_0", "Sample_1", "Sample_2"]
        test1_samples, test1_inds = _SortSamples(test1)
        test2_samples, test2_inds = _SortSamples(test2)
        test3_samples, test3_inds = _SortSamples(test3)

        assert test1_samples == ["Sample_00", "Sample_01", "Sample_02"]
        assert test1_inds == [1, 2, 0]
        assert test2_samples == ["Sample_1", "Sample_2", "Sample_3"]
        assert test2_inds == [2, 1, 0]
        assert test3_samples == ["Sample_0", "Sample_1", "Sample_2"]
        assert test3_inds == [0, 1, 2]

    def test_overlapping_samples(self):
        # Test the GetOverlappingSamples function
        snps = GenotypesVCF(None)
        strs = GenotypesTR(None)

        # Test 1 No Matching
        snps.samples = ["Sample_02", "Sample_00", "Sample_01"]
        strs.samples = ["Sample_3", "Sample_2", "Sample_1"]
        snp_inds, str_inds = GetOverlappingSamples(snps, strs)
        assert snp_inds == [] and str_inds == []

        # Test 2 All Matching
        snps.samples = ["Sample_2", "Sample_3", "Sample_1"]
        strs.samples = ["Sample_3", "Sample_2", "Sample_1"]
        snp_inds, str_inds = GetOverlappingSamples(snps, strs)
        assert snp_inds == [2, 0, 1] and str_inds == [2, 1, 0]

        # Test 3 SNPs and STRs incremented
        snps.samples = ["Sample_2", "Sample_03", "Sample_01", "Sample_4"]
        strs.samples = ["Sample_3", "Sample_2", "Sample_1", "Sample_4"]
        snp_inds, str_inds = GetOverlappingSamples(snps, strs)
        assert snp_inds == [0, 3] and str_inds == [1, 3]

        # Test 4 Uneven Sample Lists
        snps.samples = ["Sample_2", "Sample_03", "Sample_01", "Sample_4", "Sample_5"]
        strs.samples = ["Sample_3", "Sample_2", "Sample_4"]
        snp_inds, str_inds = GetOverlappingSamples(snps, strs)
        assert snp_inds == [0, 3] and str_inds == [1, 2]

    def test_gt_filter(self):
        miss = np.iinfo(np.uint8).max
        gt1 = np.array(
            [[miss - 1, 0], [1, 2], [1, miss], [miss, 1], [125, 160], [2, 2], [3, 4]]
        )
        gt2 = np.array(
            [[0, 0], [miss, miss], [4, miss - 1], [0, 0], [1, 0], [2, 1], [4, 4]]
        )
        gt1, gt2 = _FilterGts(gt1, gt2, log)
        assert np.array_equal(gt1, [285, 4, 7]) and np.array_equal(gt2, [1, 3, 8])

        gt1 = np.array([[miss, 0]])
        gt2 = np.array([[miss - 1, 1]])
        gt1, gt2 = _FilterGts(gt1, gt2, log)
        assert np.array_equal(gt1, []) and np.array_equal(gt2, [])

    def test_ld(self):
        # load expected for all tests
        expected = self._ld_expected()

        # create snp/str variants to compare
        snp1_gt = np.array(
            [[1, 0], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0], [1, 0], [1, 1], [1, 0]]
        )
        snp2_gt = np.array(
            [[1, 0], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0], [1, 0], [1, 0], [1, 1]]
        )
        str1_gt = np.array(
            [[1, 2], [1, 0], [0, 0], [0, 2], [5, 2], [0, 0], [3, 1], [4, 1], [6, 0]]
        )
        str2_gt = np.array(
            [[2, 1], [0, 1], [1, 0], [0, 0], [3, 1], [1, 1], [2, 1], [6, 0], [7, 1]]
        )

        # Calculate expected
        r1 = np.round(ComputeLD(snp1_gt, snp2_gt, "Pearson", log)[1], 4)
        r2 = np.round(ComputeLD(snp1_gt, str1_gt, "Pearson", log)[1], 4)
        r3 = np.round(ComputeLD(str1_gt, str2_gt, "Pearson", log)[1], 4)
        r4 = np.round(ComputeLD(snp1_gt, snp2_gt, "Exact", log)[1], 4)

        assert [r1, r2, r3, r4] == expected

    def test_invalid_stats_vcf(self):
        clump_p1 = 0.0001
        clump_p2 = 0.01
        clump_snp_field = "SNP"
        clump_field = "P"
        clump_chrom_field = "CHR"
        clump_pos_field = "POS"
        clump_kb = (250,)
        clump_r2 = 0.5
        LD_type = "Pearson"
        out = "NA"
        summstats_snps = "fake/path"
        summstats_strs = None
        gts_snps = None
        gts_strs = None

        def _get_error(f1, f2, ftype):
            if ftype == "SNPs":
                error = (
                    f"One of summstats-snps {f1} and gts-snps {f2} is "
                    "not present. Please ensure both have been inputted correctly."
                )
            else:
                error = (
                    f"One of summstats-strs {f1} and gts-strs {f2} is "
                    "not present. Please ensure both have been inputted correctly."
                )
            return error

        # gts_snps None
        with pytest.raises(Exception) as e:
            clumpstr(
                summstats_snps,
                summstats_strs,
                gts_snps,
                gts_strs,
                clump_p1,
                clump_p2,
                clump_snp_field,
                clump_field,
                clump_chrom_field,
                clump_pos_field,
                clump_kb,
                clump_r2,
                LD_type,
                out,
                log,
            )
        assert (str(e.value)) == _get_error(summstats_snps, gts_snps, "SNPs")

        # Summstats_snps None
        summstats_snps = None
        gts_snps = "fake/path"
        with pytest.raises(Exception) as e:
            clumpstr(
                summstats_snps,
                summstats_strs,
                gts_snps,
                gts_strs,
                clump_p1,
                clump_p2,
                clump_snp_field,
                clump_field,
                clump_chrom_field,
                clump_pos_field,
                clump_kb,
                clump_r2,
                LD_type,
                out,
                log,
            )
        assert (str(e.value)) == _get_error(summstats_snps, gts_snps, "SNPs")

        # gts_strs None
        summstats_snps = "fake/path"
        gts_snps = "fake/path"
        summstats_strs = "fake/path"
        with pytest.raises(Exception) as e:
            clumpstr(
                summstats_snps,
                summstats_strs,
                gts_snps,
                gts_strs,
                clump_p1,
                clump_p2,
                clump_snp_field,
                clump_field,
                clump_chrom_field,
                clump_pos_field,
                clump_kb,
                clump_r2,
                LD_type,
                out,
                log,
            )
        assert (str(e.value)) == _get_error(summstats_strs, gts_strs, "STRs")

        # summstats_strs None
        summstats_strs = None
        gts_strs = "fake/path"
        with pytest.raises(Exception) as e:
            clumpstr(
                summstats_snps,
                summstats_strs,
                gts_snps,
                gts_strs,
                clump_p1,
                clump_p2,
                clump_snp_field,
                clump_field,
                clump_chrom_field,
                clump_pos_field,
                clump_kb,
                clump_r2,
                LD_type,
                out,
                log,
            )
        assert (str(e.value)) == _get_error(summstats_strs, gts_strs, "STRs")


class TestClumpCLI:
    def test_snps(self, capfd):
        out_file = Path("test_snps.clump")
        cmd = (
            f"clump --summstats-snps {DATADIR}/test_snpstats.linear "
            f"--gts-snps {DATADIR}/simple.vcf "
            "--clump-id-field ID "
            "--clump-chrom-field CHROM "
            "--clump-pos-field POS "
            "--ld Exact "
            "--verbosity DEBUG "
            "--out test_snps.clump"
        )

        runner = CliRunner()
        result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
        captured = capfd.readouterr()
        assert result.exit_code == 0
        out_file.unlink()

    def test_strs(self, capfd):
        out_file = Path("test_strs.clump")
        cmd = (
            f"clump --summstats-strs {DATADIR}/test_strstats.linear "
            f"--gts-strs {DATADIR}/simple_tr.vcf "
            "--clump-id-field ID "
            "--clump-chrom-field CHROM "
            "--clump-pos-field POS "
            "--verbosity DEBUG "
            "--out test_strs.clump"
        )
        runner = CliRunner()
        result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
        assert result.exit_code == 0
        out_file.unlink()

    def test_snps_strs(self, capfd):
        out_file = Path("test_snps_strs.clump")
        cmd = (
            f"clump --summstats-snps {DATADIR}/test_snpstats.linear "
            f"--gts-snps {DATADIR}/simple.vcf "
            f"--summstats-strs {DATADIR}/test_strstats.linear "
            f"--gts-strs {DATADIR}/simple_tr.vcf "
            "--clump-id-field ID "
            "--clump-chrom-field CHROM "
            "--clump-pos-field POS "
            "--verbosity DEBUG "
            "--out test_snps_strs.clump"
        )
        runner = CliRunner()
        result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
        assert result.exit_code == 0

        out_file.unlink()
