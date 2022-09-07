import os
from pathlib import Path

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn

from haptools.transform import (
    HaplotypeAncestry,
    HaplotypesAncestry,
    GenotypesAncestry,
)
from .test_data import TestGenotypesRefAlt
from haptools.data import Variant, GenotypesRefAlt


DATADIR = Path(__file__).parent.joinpath("data")


class TestGenotypesAncestry:
    file = DATADIR.joinpath("simple-ancestry.vcf")

    def _get_expected_ancestry(self):
        # create a matrix with shape: samples x SNPs x strands
        expected = np.zeros(40).reshape((5, 4, 2)).astype(np.uint8)
        expected[0, 0, 1] = 1
        expected[1, 1, 1] = 1
        expected[2, 2, 0] = 1
        expected[3, 3, 0] = 1
        expected[4, 1, 1] = 1
        expected[3, 1, :] = 1
        expected[3, 0, 0] = 2
        expected[0, 2, 1] = 2
        expected[-1, -1, :] = 2
        return expected

    def _get_fake_genotypes(self):
        base_gts = TestGenotypesRefAlt()._get_fake_genotypes_refalt(with_phase=True)
        # copy all of the fields
        gts = GenotypesAncestry(fname=None)
        gts.data = base_gts.data
        gts.samples = base_gts.samples
        gts.variants = base_gts.variants
        # add additional ancestry data
        gts.ancestry = self._get_expected_ancestry()
        gts.ancestry_labels = {"YRI": 0, "CEU": 1, "ASW": 2}
        return gts

    def test_load_genotypes_ancestry(self, caplog):
        expected = self._get_fake_genotypes()

        # can we load the data from the VCF?
        gts = GenotypesAncestry(self.file)
        gts.read()
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.ancestry_labels == expected.ancestry_labels
        np.testing.assert_allclose(gts.ancestry, expected.ancestry)
        assert gts.samples == expected.samples

    def test_load_genotypes_iterate(self, caplog):
        expected = self._get_expected_ancestry().transpose((1, 0, 2))

        # can we load the data from the VCF?
        gts = GenotypesAncestry(self.file)
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.ancestry, expected[idx])

    def test_load_genotypes_discard_multiallelic(self):
        gts = self._get_fake_genotypes()

        # make a copy for later
        data_copy = gts.data.copy().astype(np.bool_)
        ancestry_copy = gts.ancestry.copy()
        variant_shape = list(gts.variants.shape)
        variant_shape[0] -= 1

        # force one of the SNPs to have more than one allele and check that it gets dicarded
        gts.data[1, 1, 1] = 2
        gts.check_biallelic(discard_also=True)

        data_copy_without_biallelic = np.delete(data_copy, [1], axis=1)
        ancestry_copy_without_biallelic = np.delete(ancestry_copy, [1], axis=1)
        np.testing.assert_equal(gts.data, data_copy_without_biallelic)
        np.testing.assert_equal(gts.ancestry, ancestry_copy_without_biallelic)
        assert gts.variants.shape == tuple(variant_shape)

    def test_subset_genotypes(self):
        gts = self._get_fake_genotypes()

        # subset to just the samples we want
        expected_data = gts.data[:3]
        expected_ancestry = gts.ancestry[:3]
        expected_variants = gts.variants
        samples = ("HG00096", "HG00097", "HG00099")
        gts_sub = gts.subset(samples=samples)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset to just the variants we want
        expected_data = gts.data[:, [1, 2]]
        expected_ancestry = gts.ancestry[:, [1, 2]]
        expected_variants = gts.variants[[1, 2]]
        variants = ("1:10116:A:G", "1:10117:C:A")
        gts_sub = gts.subset(variants=variants)
        assert gts_sub.samples == gts.samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset both: samples and variants
        expected_data = gts.data[[3, 4], :2]
        expected_ancestry = gts.ancestry[[3, 4], :2]
        expected_variants = gts.variants[:2]
        samples = ("HG00100", "HG00101")
        variants = ("1:10114:T:C", "1:10116:A:G")
        gts_sub = gts.subset(samples=samples, variants=variants)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

    @pytest.mark.xfail(reason="not implemented yet")
    def test_write_genotypes(self):
        assert False

    @pytest.mark.xfail(reason="not implemented yet")
    def test_write_genotypes_unphased(self):
        assert False


class TestHaplotypesAncestry:
    def _get_dummy_haps(self):
        # create three haplotypes
        haplotypes = {
            "H1": HaplotypeAncestry(
                chrom="1", start=10114, end=8, id="H1", ancestry="YRI"
            ),
            "H2": HaplotypeAncestry(
                chrom="1", start=10114, end=10119, id="H2", ancestry="YRI"
            ),
            "H3": HaplotypeAncestry(
                chrom="1", start=10116, end=10119, id="H3", ancestry="ASW"
            ),
        }
        haplotypes["H1"].variants = (
            Variant(start=10114, end=10115, id="1:10114:T:C", allele="T"),
            Variant(start=10116, end=10117, id="1:10116:A:G", allele="G"),
        )
        haplotypes["H2"].variants = (
            Variant(start=10114, end=10115, id="1:10114:T:C", allele="C"),
            Variant(start=10117, end=10118, id="1:10117:C:A", allele="C"),
        )
        haplotypes["H3"].variants = (
            Variant(start=10116, end=10117, id="1:10116:A:G", allele="A"),
            Variant(start=10117, end=10118, id="1:10117:C:A", allele="A"),
        )
        haps = HaplotypesAncestry(fname=None)
        haps.data = haplotypes
        return haps

    def test_hap_transform(self):
        expected = np.array(
            [
                [0, 0],
                [0, 0],
                [1, 1],
                [0, 0],
                [0, 0],
            ],
            dtype=np.uint8,
        )

        hap = list(self._get_dummy_haps().data.values())[0]
        gens = TestGenotypesAncestry()._get_fake_genotypes()
        gens.check_phase()
        hap_gt = hap.transform(gens)
        np.testing.assert_allclose(hap_gt, expected)

    def test_haps_transform(self):
        expected = np.array(
            [
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
                [[1, 0], [0, 1], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 1], [0, 0]],
            ],
            dtype=np.uint8,
        )

        haps = self._get_dummy_haps()
        gens = TestGenotypesAncestry()._get_fake_genotypes()
        gens.check_phase()
        gens.data[[2, 4], 0, 1] = 1
        gens.data[[1, 4], 2, 0] = 1
        hap_gt = GenotypesRefAlt(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)
        return hap_gt
