import os
from pathlib import Path

import pytest
import numpy as np
from cyvcf2 import VCF

from haptools.data.genotypes import GenotypesRefAlt
from haptools.data.genotypes import GenotypesTR
from haptools.clump import (
    LoadVariant,
    Variant
)

DATADIR = Path(__file__).parent.joinpath("data")

def test_loading_snps():
    gts_snps = DATADIR.joinpath("outvcf_test.vcf.gz")
    snpgts = GenotypesRefAlt.load(str(gts_snps))
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
        vargts = LoadVariant(var, snpgts, strgts)
        assert len(vargts) == snpgts.data.shape[0]
        assert np.array_equal(vargts, answer)
    return

def test_loading_strs():
    return

def test_loading_snps_strs():
    # TODO test overlap functionality and that correct variant is found
    return