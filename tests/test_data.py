import os
from pathlib import Path
from dataclasses import dataclass, field

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn

from haptools.sim_phenotype import Haplotype as HaptoolsHaplotype
from haptools.data import (
    Extra,
    Variant,
    Haplotype,
    Genotypes,
    Phenotypes,
    Covariates,
    Haplotypes,
    GenotypesRefAlt,
    GenotypesPLINK,
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

    def _get_fake_genotypes(self):
        gts = Genotypes(fname=None)
        gts.data = self._get_expected_genotypes()
        gts.variants = np.array(
            [
                ("1:10114:T:C", "1", 10114, 0),
                ("1:10116:A:G", "1", 10116, 0.6),
                ("1:10117:C:A", "1", 10117, 0),
                ("1:10122:A:G", "1", 10122, 0),
            ],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint32),
                ("aaf", np.float64),
            ],
        )
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.check_phase()
        return gts

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

        # force one of the samples to have a missing GT and check that we get an error
        gts.data[1, 1, 1] = -1
        with pytest.raises(ValueError) as info:
            gts.check_missing()
        assert (
            str(info.value)
            == "Genotype with ID 1:10116:A:G at POS 1:10116 is missing for sample"
            " HG00097"
        )
        gts.data[1, 1, 1] = 1

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
            == "Variant with ID 1:10116:A:G at POS 1:10116 is unphased for sample"
            " HG00097"
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

    def test_load_genotypes_example(self):
        samples = (
            "HG00096",
            "HG00097",
            "HG00099",
            "HG00100",
            "HG00101",
            "HG00102",
            "HG00103",
            "HG00105",
            "HG00106",
            "HG00107",
            "HG00108",
            "HG00109",
            "HG00110",
            "HG00111",
            "HG00112",
            "HG00113",
            "HG00114",
            "HG00115",
            "HG00116",
            "HG00117",
            "HG00118",
            "HG00119",
            "HG00120",
            "HG00121",
            "HG00122",
        )
        gts = Genotypes(DATADIR.joinpath("example.vcf.gz"))
        gts.read()
        assert gts.samples[:25] == samples

    def test_load_genotypes_iterate(self, caplog):
        expected = self._get_expected_genotypes().transpose((1, 0, 2))
        samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR.joinpath("simple.vcf"))
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data, expected[idx])
        assert gts.samples == samples

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

        # subset to just one of the variants
        expected = expected[:, [1]]

        gts = Genotypes(DATADIR.joinpath("simple.vcf.gz"))
        samples = ["HG00097", "HG00100"]
        variants = {"1:10117:C:A"}
        gts.read(region="1:10115-10117", samples=samples, variants=variants)
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == tuple(samples)

    def test_subset_genotypes(self):
        gts = self._get_fake_genotypes()

        # subset to just the samples we want
        expected_data = gts.data[:3]
        expected_variants = gts.variants
        samples = ("HG00096", "HG00097", "HG00099")
        gts_sub = gts.subset(samples=samples)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset to just the variants we want
        expected_data = gts.data[:, [1, 2]]
        expected_variants = gts.variants[[1, 2]]
        variants = ("1:10116:A:G", "1:10117:C:A")
        gts_sub = gts.subset(variants=variants)
        assert gts_sub.samples == gts.samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset both: samples and variants
        expected_data = gts.data[[3, 4], :2]
        expected_variants = gts.variants[:2]
        samples = ("HG00100", "HG00101")
        variants = ("1:10114:T:C", "1:10116:A:G")
        gts_sub = gts.subset(samples=samples, variants=variants)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)


class TestGenotypesPLINK:
    def _get_fake_genotypes_plink(self):
        pgenlib = pytest.importorskip("pgenlib")
        gts_ref_alt = TestGenotypesRefAlt()._get_fake_genotypes_refalt()
        gts = GenotypesPLINK(gts_ref_alt.fname)
        gts.data = gts_ref_alt.data
        gts.samples = gts_ref_alt.samples
        gts.variants = gts_ref_alt.variants
        return gts

    def test_load_genotypes(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR.joinpath("simple.pgen"))
        gts.read()
        gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "ref", "alt"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_iterate(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR.joinpath("simple.pgen"))

        # check that everything matches what we expected
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :2], expected.data[:, idx])
            for col in ("chrom", "pos", "id", "ref", "alt"):
                assert line.variants[col] == expected.variants[col][idx]
        assert gts.samples == expected.samples

    def test_load_genotypes_subset(self):
        expected = self._get_fake_genotypes_plink()

        # subset for the region we want
        expected_data = expected.data[:, 1:3]

        # can we load the data from the VCF?
        gts = GenotypesPLINK(DATADIR.joinpath("simple.pgen"))
        gts.read(region="1:10115-10117")
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == expected.samples

        # subset for just the samples we want
        expected_data = expected_data[[1, 3]]

        gts = GenotypesPLINK(DATADIR.joinpath("simple.pgen"))
        samples = [expected.samples[1], expected.samples[3]]
        gts.read(region="1:10115-10117", samples=samples)
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == tuple(samples)

        # subset to just one of the variants
        expected_data = expected_data[:, [1]]

        gts = GenotypesPLINK(DATADIR.joinpath("simple.pgen"))
        variants = {"1:10117:C:A"}
        gts.read(region="1:10115-10117", samples=samples, variants=variants)
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == tuple(samples)

    def test_write_genotypes(self):
        gts = self._get_fake_genotypes_plink()

        fname = DATADIR.joinpath("test_write.pgen")
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts.read()
        new_gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "ref", "alt"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_genotypes_unphased(self):
        gts = self._get_fake_genotypes_plink()
        # add phasing information back
        gts.data = np.dstack((gts.data, np.ones(gts.data.shape[:2], dtype=np.uint8)))
        gts.data[:2, 1, 2] = 0

        fname = DATADIR.joinpath("test_unphased.pgen")
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts.read()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "ref", "alt"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()


class TestPhenotypes:
    def _get_expected_phenotypes(self):
        # create a phen matrix with shape: samples x phenotypes
        expected = np.array(
            [
                [1, 1, 2, 2, 0],
                [3, 6, 8, 1, 4],
            ],
            dtype="float64",
        ).T
        return expected

    def _get_fake_phenotypes(self):
        gts = Phenotypes(fname=None)
        gts.data = self._get_expected_phenotypes()
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.names = ("height", "bmi")
        return gts

    def test_load_phenotypes(self, caplog):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR.joinpath("simple.pheno"))
        phens.read()
        np.testing.assert_allclose(phens.data, expected)
        assert phens.samples == expected_phen.samples

        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        phens.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

        expected = (expected - np.mean(expected, axis=0)) / np.std(expected, axis=0)
        phens.standardize()
        np.testing.assert_allclose(phens.data, expected)

    def test_load_phenotypes_iterate(self):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data
        samples = expected_phen.samples

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR.joinpath("simple.pheno"))
        for idx, line in enumerate(phens):
            np.testing.assert_allclose(line.data, expected[idx])
            assert line.samples == samples[idx]

    def test_load_phenotypes_subset(self):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data
        samples = ["HG00097", "HG00100"]

        # subset for just the samples we want
        expected = expected[[1, 3]]

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR.joinpath("simple.pheno"))
        phens.read(samples=samples)
        np.testing.assert_allclose(phens.data, expected)
        assert phens.samples == tuple(samples)

    def test_write_phenotypes(self):
        expected_phen = self._get_fake_phenotypes()

        # first, we write the data
        expected_phen.fname = DATADIR.joinpath("test.pheno")
        expected_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(expected_phen.fname)
        result.read()
        np.testing.assert_allclose(expected_phen.data, result.data)
        assert expected_phen.names == result.names
        assert expected_phen.samples == result.samples

        # let's clean up after ourselves and delete the file
        os.remove(str(expected_phen.fname))

    def test_append_phenotype(self):
        expected1 = self._get_fake_phenotypes()
        expected2 = self._get_fake_phenotypes()

        # add a phenotype called "test" to the set of phenotypes
        new_phen = np.array([5, 2, 8, 0, 3], dtype=expected2.data.dtype)
        expected2.append("test", new_phen)

        # did it get added correctly?
        assert expected1.data.shape[1] == expected2.data.shape[1] - 1
        assert len(expected1.names) == len(expected2.names) - 1
        assert expected2.names[2] == "test"
        assert len(expected1.samples) == len(expected2.samples)
        np.testing.assert_allclose(expected1.data, expected2.data[:, :2])
        np.testing.assert_allclose(expected2.data[:, 2], new_phen)


class TestCovariates:
    def _get_expected_covariates(self):
        # create a covar matrix with shape: samples x covariates
        expected = np.array([(0, 4), (1, 20), (1, 33), (0, 15), (0, 78)])
        return expected

    def _get_fake_covariates(self):
        gts = Phenotypes(fname=None)
        gts.data = self._get_expected_covariates()
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.names = ("sex", "age")
        return gts

    def test_load_covariates(self, caplog):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = expected_covar.samples

        # can we load the data from the covariates file?
        covars = Covariates(DATADIR.joinpath("simple.covar"))
        covars.read()
        np.testing.assert_allclose(covars.data, expected)
        assert covars.samples == samples
        assert covars.names == ("sex", "age")

        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        covars.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

    def test_load_covariates_iterate(self):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = expected_covar.samples

        # can we load the data from the covariates file?
        covars = Covariates(DATADIR.joinpath("simple.covar"))
        for idx, line in enumerate(covars):
            np.testing.assert_allclose(line.data, expected[idx])
            assert line.samples == samples[idx]
        assert covars.names == ("sex", "age")

    def test_load_covariates_subset(self):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = ["HG00097", "HG00100"]

        # subset for just the samples we want
        expected = expected[[1, 3]]

        # can we load the data from the covariate file?
        covars = Covariates(DATADIR.joinpath("simple.covar"))
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

    def _get_dummy_haps(self):
        # create three haplotypes
        haplotypes = {
            "H1": Haplotype(chrom="1", start=10114, end=8, id="H1"),
            "H2": Haplotype(chrom="1", start=10114, end=10119, id="H2"),
            "H3": Haplotype(chrom="1", start=10116, end=10119, id="H3"),
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
        haps = Haplotypes(fname=None)
        haps.data = haplotypes
        return haps

    def test_load(self):
        expected = self._basic_haps()

        # can we load this data from the hap file?
        haps = Haplotypes.load(DATADIR.joinpath("basic.hap"))
        assert expected == haps.data

        # also check the indexed file
        # it should be the same
        haps = Haplotypes.load(DATADIR.joinpath("basic.hap.gz"))
        assert expected == haps.data

    def test_read_subset(self):
        expected = {}
        expected["chr21.q.3365*1"] = self._basic_haps()["chr21.q.3365*1"]

        haps = Haplotypes(DATADIR.joinpath("basic.hap"))
        # this shouldn't fail anymore as of version 0.1.0
        haps.read(haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

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

    def _get_writable_haplotypes(self):
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
        return expected

    def test_write_extras(self):
        # what do we expect to see from the test.hap file?
        expected = self._get_writable_haplotypes()

        haps = Haplotypes(DATADIR.joinpath("test.hap"), HaptoolsHaplotype)
        haps.data = expected
        haps.write()

        haps.data = None
        haps.read()
        assert expected == haps.data

        # remove the file
        os.remove("tests/data/test.hap")

    def test_write_plus_extra(self):
        @dataclass
        class HaplotypePlusExtra(HaptoolsHaplotype):
            """
            A haplotype with an additional, unnecessary extra field

            Properties and functions are shared with the HaptoolsHaplotype object
            """

            score: float
            _extras: tuple = field(
                repr=False,
                init=False,
                default=(
                    Extra("ancestry", "s", "Local ancestry"),
                    Extra("score", ".2f", "Score for a thing"),
                    Extra("beta", ".2f", "Effect size in linear model"),
                ),
            )
        # what do we want to write to the test.hap file?
        expected = {
            "chr21.q.3365*1": HaplotypePlusExtra(
                "21", 26928472, 26941960, "chr21.q.3365*1", "ASW", 0.73, 0.40
            ),
            "chr21.q.3365*10": HaplotypePlusExtra(
                "21", 26938989, 26941960, "chr21.q.3365*10", "CEU", 0.30, 0.28
            ),
            "chr21.q.3365*11": HaplotypePlusExtra(
                "21", 26938353, 26938989, "chr21.q.3365*11", "MXL", 0.49, 0.84
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            expected[hap_id].variants = hap.variants

        haps = Haplotypes(DATADIR.joinpath("test.hap"), HaplotypePlusExtra)
        haps.data = expected
        haps.write()

        haps = Haplotypes(DATADIR.joinpath("test.hap"), HaptoolsHaplotype)
        haps.read()
        assert haps.data == self._get_writable_haplotypes()

        # remove the file
        os.remove("tests/data/test.hap")

    def test_hap_transform(self):
        expected = np.array(
            [
                [0, 1],
                [0, 1],
                [1, 1],
                [1, 1],
                [0, 0],
            ],
            dtype=np.uint8,
        )

        hap = list(self._get_dummy_haps().data.values())[0]
        gens = TestGenotypesRefAlt()._get_fake_genotypes_refalt()
        hap_gt = hap.transform(gens)
        np.testing.assert_allclose(hap_gt, expected)

    def test_haps_transform(self):
        expected = np.array(
            [
                [[0, 1], [0, 0], [0, 0]],
                [[0, 1], [0, 0], [1, 0]],
                [[1, 0], [0, 1], [0, 0]],
                [[1, 1], [0, 0], [0, 0]],
                [[0, 0], [0, 1], [1, 0]],
            ],
            dtype=np.uint8,
        )

        haps = self._get_dummy_haps()
        gens = TestGenotypesRefAlt()._get_fake_genotypes_refalt()
        gens.data[[2, 4], 0, 1] = 1
        gens.data[[1, 4], 2, 0] = 1
        hap_gt = GenotypesRefAlt(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)
        return hap_gt

    def test_hap_gt_write(self):
        fname = DATADIR.joinpath("simple_haps.vcf")

        hap_gt = self.test_haps_transform()
        hap_gt.fname = fname
        expected_data = hap_gt.data
        expected_samples = hap_gt.samples
        hap_gt.write()

        hap_gt.data = None
        hap_gt.read()
        hap_gt.check_phase()
        np.testing.assert_allclose(hap_gt.data, expected_data)
        assert hap_gt.samples == expected_samples
        assert len(hap_gt.variants) == hap_gt.data.shape[1]

        # remove the file
        os.remove(str(fname))


class TestGenotypesRefAlt:
    def _get_fake_genotypes_refalt(self):
        base_gts = TestGenotypes()._get_fake_genotypes()
        # copy all of the fields
        gts = GenotypesRefAlt(fname=None)
        gts.data = base_gts.data
        gts.samples = base_gts.samples
        # add additional ref and alt alleles
        ref_alt = np.array(
            [
                ("T", "C"),
                ("A", "G"),
                ("C", "A"),
                ("A", "G"),
            ],
            dtype=[
                ("ref", "U100"),
                ("alt", "U100"),
            ],
        )
        # see https://stackoverflow.com/a/5356137
        gts.variants = rfn.merge_arrays((base_gts.variants, ref_alt), flatten=True)
        return gts

    def test_read_ref_alt(self):
        # simple.vcf
        gts_ref_alt_read = GenotypesRefAlt(DATADIR.joinpath("simple.vcf"))
        gts_ref_alt_read.read()
        expected = np.array(
            [
                ("T", "C"),
                ("A", "G"),
                ("C", "A"),
                ("A", "G"),
            ],
            dtype=gts_ref_alt_read.variants[["ref", "alt"]].dtype,
        )
        for i, x in enumerate(expected):
            assert gts_ref_alt_read.variants[["ref", "alt"]][i] == x

        # example.vcf.gz
        gts_ref_alt = GenotypesRefAlt(DATADIR.joinpath("example.vcf.gz"))
        gts_ref_alt.read()
        expected = np.array(
            [
                ("C", "A"),
                ("T", "C"),
                ("G", "A"),
                ("T", "C"),
                ("A", "G"),
                ("A", "G"),
                ("T", "G"),
                ("T", "A"),
                ("G", "A"),
                ("T", "C"),
                ("C", "G"),
                ("A", "G"),
                ("T", "C"),
                ("T", "C"),
                ("C", "T"),
            ],
            dtype=gts_ref_alt.variants[["ref", "alt"]].dtype,
        )
        for i, x in enumerate(expected):
            assert gts_ref_alt.variants[["ref", "alt"]][i] == x

    def test_write_ref_alt(self):
        # strategy is to read in the file, write it, and then read again
        # read genotypes
        gts_ref_alt_write = GenotypesRefAlt(DATADIR.joinpath("simple.vcf"))
        gts_ref_alt_write.read()
        gts_ref_alt_write.check_phase()
        # write file to new file
        fname = DATADIR.joinpath("test.vcf")
        gts_ref_alt_write.fname = fname
        gts_ref_alt_write.write()
        # read again
        gts_ref_alt_write.read()
        gts_ref_alt_write.check_phase()
        # compare samples, data, variants, and ref/alt for equality
        expected = self._get_fake_genotypes_refalt()
        assert gts_ref_alt_write.samples == expected.samples
        np.testing.assert_allclose(gts_ref_alt_write.data, expected.data)

        for i, x in enumerate(expected.variants):
            # dtype.names gives us names of columns in the array
            for col in expected.variants.dtype.names:
                # each row is a variant
                # index into col, i gets specific variant, x iterates thru
                assert gts_ref_alt_write.variants[col][i] == x[col]

        os.remove(str(fname))
