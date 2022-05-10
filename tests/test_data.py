import os
from pathlib import Path

import pytest
import numpy as np

from haptools.haplotype import HaptoolsHaplotype
from haptools.data import (
    Genotypes,
    Phenotypes,
    Covariates,
    Haplotypes,
    Variant,
    Haplotype,
)


DATADIR = Path(__file__).parent.joinpath("data")


class TestGenotypes:
    def _get_expected_genotypes(self):
        # create a GT matrix with shape: samples x SNPs x (strands+phase)
        expected = np.zeros(60).reshape((5, 4, 3)).astype(np.uint8)
        expected[:4, 1, 1] = 1
        expected[2:4, 1, 0] = 1
        expected[:, :, 2] = 1
        return expected


    def test_load_genotypes(self, caplog):
        expected = self._get_expected_genotypes()

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR.joinpath("simple.vcf"))
        gts.read()
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")


        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        gts.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

        # force one of the SNPs to have more than one allele and check that we get an error
        gts.data[1, 1, 1] = 2
        with pytest.raises(ValueError) as info:
            gts.check_biallelic()
        assert (
            str(info.value)
            == "Variant with ID 1:10116:A:G at POS 1:10116 is multiallelic for sample"
            " HG00097"
        )
        gts.data[1, 1, 1] = 1

        # check biallelic-ness and convert to bool_
        gts.check_biallelic()
        expected = expected.astype(np.bool_)
        np.testing.assert_allclose(gts.data, expected)

        # force one of the het SNPs to be unphased and check that we get an error message
        gts.data[1, 1, 2] = 0
        with pytest.raises(ValueError) as info:
            gts.check_phase()
        assert (
            str(info.value)
            == "Variant with ID 1:10116:A:G at POS 1:10116 is unphased for sample HG00097"
        )
        gts.data[1, 1, 2] = 1

        # check phase and remove the phase axis
        gts.check_phase()
        expected = expected[:, :, :2]
        np.testing.assert_allclose(gts.data, expected)

        # try to check phase again - it should warn b/c we've already done it before
        caplog.clear()
        gts.check_phase()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

        # convert the matrix of alt allele counts to a matrix of minor allele counts
        assert gts.variants["aaf"][1] == 0.6
        gts.to_MAC()
        expected[:, 1, :] = ~expected[:, 1, :]
        np.testing.assert_allclose(gts.data, expected)
        assert gts.variants["maf"][1] == 0.4

        # try to do the MAC conversion again - it should warn b/c we've already done it
        caplog.clear()
        gts.to_MAC()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"


    def test_load_genotypes_iterate(self, caplog):
        expected = self._get_expected_genotypes().transpose((1, 0, 2))
        samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR.joinpath("simple.vcf"))
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data, expected[idx])
            assert line.samples == samples


    def test_load_genotypes_discard_multiallelic(self):
        expected = self._get_expected_genotypes()

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR.joinpath("simple.vcf"))
        gts.read()

        # make a copy for later
        data_copy = gts.data.copy().astype(np.bool_)
        variant_shape = list(gts.variants.shape)
        variant_shape[0] -= 1

        # force one of the SNPs to have more than one allele and check that it gets dicarded
        gts.data[1, 1, 1] = 2
        gts.check_biallelic(discard_also=True)

        data_copy_without_biallelic = np.delete(data_copy, [1], axis=1)
        np.testing.assert_equal(gts.data, data_copy_without_biallelic)
        assert gts.variants.shape == tuple(variant_shape)


    def test_load_genotypes_subset(self):
        expected = self._get_expected_genotypes()

        # subset for the region we want
        expected = expected[:, 1:3]

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR.joinpath("simple.vcf.gz"))
        gts.read(region="1:10115-10117")
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        # subset for just the samples we want
        expected = expected[[1, 3]]

        gts = Genotypes(DATADIR.joinpath("simple.vcf.gz"))
        samples = ["HG00097", "HG00100"]
        gts.read(region="1:10115-10117", samples=samples)
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == tuple(samples)


def test_load_phenotypes(caplog):
    # create a phenotype vector with shape: num_samples x 1
    expected = np.array([1, 1, 2, 2, 0])

    # can we load the data from the phenotype file?
    phens = Phenotypes(DATADIR.joinpath("simple.tsv"))
    phens.read()
    np.testing.assert_allclose(phens.data, expected)
    assert phens.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # try loading the data again - it should warn b/c we've already done it
    phens.read()
    assert len(caplog.records) == 1 and caplog.records[0].levelname == "WARNING"

    expected = (expected - np.mean(expected)) / np.std(expected)
    phens.standardize()
    np.testing.assert_allclose(phens.data, expected)


def test_load_phenotypes_iterate():
    # create a phenotype vector with shape: num_samples x 1
    expected = np.array([1, 1, 2, 2, 0])
    samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # can we load the data from the phenotype file?
    phens = Phenotypes(DATADIR.joinpath("simple.tsv"))
    for idx, line in enumerate(phens):
        np.testing.assert_allclose(line.data, expected[idx])
        assert line.samples == samples[idx]


def test_load_phenotypes_subset():
    # create a phenotype vector with shape: num_samples x 1
    expected = np.array([1, 1, 2, 2, 0])

    # subset for just the samples we want
    expected = expected[[1, 3]]

    # can we load the data from the phenotype file?
    phens = Phenotypes(DATADIR.joinpath("simple.tsv"))
    samples = ["HG00097", "HG00100"]
    phens.read(samples=samples)
    np.testing.assert_allclose(phens.data, expected)
    assert phens.samples == tuple(samples)


def test_load_covariates(caplog):
    # create a covariate vector with shape: num_samples x num_covars
    expected = np.array([(0, 4), (1, 20), (1, 33), (0, 15), (0, 78)])

    # can we load the data from the covariates file?
    covars = Covariates(DATADIR.joinpath("covars.tsv"))
    covars.read()
    np.testing.assert_allclose(covars.data, expected)
    assert covars.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
    assert covars.names == ("sex", "age")

    # try loading the data again - it should warn b/c we've already done it
    covars.read()
    assert len(caplog.records) == 1 and caplog.records[0].levelname == "WARNING"


def test_load_covariates_iterate():
    # create a covariate vector with shape: num_samples x num_covars
    expected = np.array([(0, 4), (1, 20), (1, 33), (0, 15), (0, 78)])
    samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # can we load the data from the covariates file?
    covars = Covariates(DATADIR.joinpath("covars.tsv"))
    for idx, line in enumerate(covars):
        np.testing.assert_allclose(line.data, expected[idx])
        assert line.samples == samples[idx]
        assert line.names == ("sex", "age")


def test_load_covariates_subset():
    # create a covariate vector with shape: num_samples x num_covars
    expected = np.array([(0, 4), (1, 20), (1, 33), (0, 15), (0, 78)])

    # subset for just the samples we want
    expected = expected[[1, 3]]

    # can we load the data from the covariate file?
    covars = Covariates(DATADIR.joinpath("covars.tsv"))
    samples = ["HG00097", "HG00100"]
    covars.read(samples=samples)
    np.testing.assert_allclose(covars.data, expected)
    assert covars.samples == tuple(samples)


class TestHaplotypes:
    def _basic_haps(self):
        # what do we expect to see from the basic.hap file?
        expected = {
            "chr21.q.3365*1": Haplotype("21", 26928472, 26941960, "chr21.q.3365*1"),
            "chr21.q.3365*10": Haplotype("21", 26938989, 26941960, "chr21.q.3365*10"),
            "chr21.q.3365*11": Haplotype("21", 26938353, 26938989, "chr21.q.3365*11"),
        }
        expected["chr21.q.3365*1"].variants = (
            Variant(26928472, 26928472, "21_26928472_C_A", "C"),
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26940815, 26940815, "21_26940815_T_C", "C"),
            Variant(26941960, 26941960, "21_26941960_A_G", "G"),
        )
        expected["chr21.q.3365*10"].variants = (
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
            Variant(26940815, 26940815, "21_26940815_T_C", "T"),
            Variant(26941960, 26941960, "21_26941960_A_G", "A"),
        )
        expected["chr21.q.3365*11"].variants = (
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
        )
        return expected

    def test_load(self):
        # can we load this data from the hap file?
        haps = Haplotypes.load(DATADIR.joinpath("basic.hap"))
        assert self._basic_haps() == haps.data

    def test_read_subset(self):
        expected = {}
        expected["chr21.q.3365*1"] = self._basic_haps()["chr21.q.3365*1"]

        haps = Haplotypes(DATADIR.joinpath("basic.hap"))
        # this should fail with an informative error b/c the file isn't indexed
        with pytest.raises(OSError) as info:
            haps.read(haplotypes={"chr21.q.3365*1"})
        assert len(str(info.value))

        haps = Haplotypes(DATADIR.joinpath("basic.hap.gz"))
        haps.read(haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

        haps = Haplotypes(DATADIR.joinpath("basic.hap.gz"))
        haps.read(region="21:26928472-26941960", haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

        expected = self._basic_haps()

        haps = Haplotypes(DATADIR.joinpath("basic.hap.gz"))
        haps.read(region="21:26928472-26941960")
        assert expected == haps.data

    def test_read_extras(self):
        # what do we expect to see from the simphenotype.hap file?
        expected = {
            "chr21.q.3365*1": HaptoolsHaplotype(
                "21", 26928472, 26941960, "chr21.q.3365*1", "ASW", 0.73
            ),
            "chr21.q.3365*10": HaptoolsHaplotype(
                "21", 26938989, 26941960, "chr21.q.3365*10", "CEU", 0.30
            ),
            "chr21.q.3365*11": HaptoolsHaplotype(
                "21", 26938353, 26938989, "chr21.q.3365*11", "MXL", 0.49
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            expected[hap_id].variants = hap.variants

        # can we load this data from the hap file?
        haps = Haplotypes(DATADIR.joinpath("simphenotype.hap"), HaptoolsHaplotype)
        haps.read()
        assert expected == haps.data

    def test_read_extras_large(self):
        """
        try reading a large-ish file
        """
        haps = Haplotypes(DATADIR.joinpath("example.hap.gz"), HaptoolsHaplotype)
        haps.read()
        assert len(self._basic_haps().keys() & haps.data.keys())

    def test_write(self):
        haps = Haplotypes(DATADIR.joinpath("test.hap"))
        haps.data = self._basic_haps()
        haps.write()

        haps.data = None
        haps.read()
        assert self._basic_haps() == haps.data

        # remove the file
        os.remove("tests/data/test.hap")

    def test_write_extras(self):
        # what do we expect to see from the simphenotype.hap file?
        expected = {
            "chr21.q.3365*1": HaptoolsHaplotype(
                "21", 26928472, 26941960, "chr21.q.3365*1", "ASW", 0.73
            ),
            "chr21.q.3365*10": HaptoolsHaplotype(
                "21", 26938989, 26941960, "chr21.q.3365*10", "CEU", 0.30
            ),
            "chr21.q.3365*11": HaptoolsHaplotype(
                "21", 26938353, 26938989, "chr21.q.3365*11", "MXL", 0.49
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            expected[hap_id].variants = hap.variants

        haps = Haplotypes(DATADIR.joinpath("test.hap"), HaptoolsHaplotype)
        haps.data = expected
        haps.write()

        haps.data = None
        haps.read()
        assert expected == haps.data

        # remove the file
        os.remove("tests/data/test.hap")
