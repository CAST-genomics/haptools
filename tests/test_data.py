import pytest
import numpy as np

from pathlib import Path

from haptools.data import Genotypes, Phenotypes, Covariates


DATADIR = Path(__file__).parent.joinpath("data")


def get_expected_genotypes():
    # create a GT matrix with shape: samples x SNPs x (strands+phase)
    expected = np.zeros(60).reshape((5, 4, 3)).astype(np.uint8)
    expected[:4, 1, 1] = 1
    expected[2:4, 1, 0] = 1
    expected[:, :, 2] = 1
    return expected


def test_load_genotypes():
    expected = get_expected_genotypes()

    # can we load the data from the VCF?
    gts = Genotypes(DATADIR.joinpath("simple.vcf"))
    gts.read()
    np.testing.assert_allclose(gts.data, expected)
    assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # try loading the data again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        gts.read()

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

    # try to check phase again - it should fail b/c we've already done it before
    with pytest.raises(AssertionError):
        gts.check_phase()

    # convert the matrix of alt allele counts to a matrix of minor allele counts
    assert gts.variants["aaf"][1] == 0.6
    gts.to_MAC()
    expected[:, 1, :] = ~expected[:, 1, :]
    np.testing.assert_allclose(gts.data, expected)
    assert gts.variants["maf"][1] == 0.4

    # try to do the MAC conversion again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        gts.to_MAC()


def test_load_genotypes_subset():
    expected = get_expected_genotypes()

    # subset for the region we want
    expected = expected[:, 1:3]

    # can we load the data from the VCF?
    gts = Genotypes(DATADIR.joinpath("simple.vcf.gz"))
    gts.read(region="1:10115-10117")
    np.testing.assert_allclose(gts.data, expected)
    assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # subset for just the samples we want
    expected = expected[[1, 3]]

    # can we load the data from the VCF?
    gts = Genotypes(DATADIR.joinpath("simple.vcf.gz"))
    samples = ["HG00097", "HG00100"]
    gts.read(region="1:10115-10117", samples=samples)
    np.testing.assert_allclose(gts.data, expected)
    assert gts.samples == tuple(samples)


def test_load_phenotypes():
    # create a phenotype vector with shape: num_samples x 1
    expected = np.array([1, 1, 2, 2, 0])

    # can we load the data from the phenotype file?
    phens = Phenotypes(DATADIR.joinpath("simple.tsv"))
    phens.read()
    np.testing.assert_allclose(phens.data, expected)
    assert phens.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # try loading the data again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        phens.read()

    expected = (expected - np.mean(expected)) / np.std(expected)
    phens.standardize()
    np.testing.assert_allclose(phens.data, expected)


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

def test_load_covariates():
    # create a covariate vector with shape: num_samples x num_covars
    expected = np.array([
        (0, 4),
        (1, 20),
        (1, 33),
        (0, 15),
        (0, 78)
    ])

    # can we load the data from the covariates file?
    covars = Covariates(DATADIR.joinpath("covars.tsv"))
    covars.read()
    np.testing.assert_allclose(covars.data, expected)
    assert covars.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
    assert covars.names == ('sex', 'age')

    # try loading the data again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        covars.read()

def test_load_covariates_subset():
    # create a covriate vector with shape: num_samples x num_covars
    expected = np.array([
        (0, 4),
        (1, 20),
        (1, 33),
        (0, 15),
        (0, 78)
    ])

    # subset for just the samples we want
    expected = expected[[1, 3]]

    # can we load the data from the covariate file?
    covars = Covariates(DATADIR.joinpath("covars.tsv"))
    samples = ["HG00097", "HG00100"]
    covars.read(samples=samples)
    np.testing.assert_allclose(covars.data, expected)
    assert covars.samples == tuple(samples)
