import os
from pathlib import Path

import pytest
import numpy as np
from cyvcf2 import VCF

from haptools.logging import getLogger
from haptools.data.genotypes import GenotypesVCF
from haptools.data.genotypes import GenotypesTR
from haptools.clump import (
    GetOverlappingSamples,
    _SortSamples,
    LoadVariant,
    _FilterGts,
    Variant
)

DATADIR = Path(__file__).parent.joinpath("data")
log = getLogger(name="test")

def test_loading_snps():
    gts_snps = DATADIR.joinpath("outvcf_test.vcf.gz")
    snpgts = GenotypesVCF.load(str(gts_snps))
    snpvars = [
                Variant("test1", "1", "10114", "0.05", "snp"),
                Variant("test2", "1", "59423090", "0.05", "snp"),
                Variant("test3", "2", "10122", "0.05", "snp")
              ]
    strgts = None

    answers = [
                np.array([2, 2, 0, 0, 0]),
                np.array([0, 0, 2, 2, 2]),
                np.array([0, 0, 2, 2, 2])
              ]

    for var, answer in zip(snpvars, answers):
        vargts = LoadVariant(var, snpgts, strgts, log)
        assert len(vargts) == snpgts.data.shape[0]
        assert np.array_equal(vargts, answer)


def test_sample_sorting():
    test1 = ["Sample_02", "Sample_00", "Sample_01"]
    test2 = ["Sample_3", "Sample_2", "Sample_1"]
    test3 = ["Sample_0", "Sample_1", "Sample_2"]
    test1_samples, test1_inds = _SortSamples(test1)
    test2_samples, test2_inds = _SortSamples(test2)
    test3_samples, test3_inds = _SortSamples(test3)

    assert test1_samples == ["Sample_00", "Sample_01", "Sample_02"]
    assert test1_inds == [1,2,0]
    assert test2_samples == ["Sample_1", "Sample_2", "Sample_3"]
    assert test2_inds == [2,1,0]
    assert test3_samples == ["Sample_0", "Sample_1", "Sample_2"]
    assert test3_inds == [0,1,2]


def test_overlapping_samples():
    # Test the GetOverlappingSamples function
    snps = GenotypesVCF(fname="NA")
    strs = GenotypesTR(fname="NA")

    # Test 1 No Matching
    snps.samples = ["Sample_02", "Sample_00", "Sample_01"]
    strs.samples = ["Sample_3", "Sample_2", "Sample_1"]
    snp_inds, str_inds = GetOverlappingSamples(snps, strs)
    assert snp_inds == [] and str_inds == []

    # Test 2 All Matching
    snps.samples = ["Sample_2", "Sample_3", "Sample_1"]
    strs.samples = ["Sample_3", "Sample_2", "Sample_1"]
    snp_inds, str_inds = GetOverlappingSamples(snps, strs)
    assert snp_inds == [2,0,1] and str_inds == [2,1,0]

    # Test 3 SNPs and STRs incremented
    snps.samples = ["Sample_2", "Sample_03", "Sample_01", "Sample_4"]
    strs.samples = ["Sample_3", "Sample_2", "Sample_1", "Sample_4"]
    snp_inds, str_inds = GetOverlappingSamples(snps, strs)
    assert snp_inds == [0,3] and str_inds == [1,3]

    # Test 4 Uneven Sample Lists
    snps.samples = ["Sample_2", "Sample_03", "Sample_01", "Sample_4", "Sample_5"]
    strs.samples = ["Sample_3", "Sample_2", "Sample_4"]
    snp_inds, str_inds = GetOverlappingSamples(snps, strs)
    assert snp_inds == [0,3] and str_inds == [1,2]

def test_gt_filter():
    gt1 = np.array([255, 1, 255, 3, 4, 7])
    gt2 = np.array([0, 255, 255, 1, 3, 8])
    gt1, gt2 = _FilterGts(gt1, gt2)
    assert np.array_equal(gt1, [3,4,7]) and np.array_equal(gt2, [1,3,8])

    gt1 = np.array([255])
    gt2 = np.array([255])
    gt1, gt2 = _FilterGts(gt1, gt2)
    assert np.array_equal(gt1, []) and np.array_equal(gt2, [])
