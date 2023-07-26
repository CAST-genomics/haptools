import os
from pathlib import Path

import pytest

from . import test_data
from haptools import val_hapfile
from haptools import data

DATADIR = Path(__file__).parent.joinpath("data").joinpath("hapfiles")


def _generate_fake_haps():
    haps_ = test_data.TestHaplotypes()
    haps = haps_._get_dummy_haps()
    haps.fname = Path(DATADIR / "valhap_test_data.hap")
    haps.write()


def _generate_fake_vars():
    vars_ = test_data.TestGenotypesPLINK()
    vars = vars_._get_fake_genotypes_plink()
    vars.fname = Path(DATADIR / "valhap_test_data.plink")
    vars.write_variants()


def test_generated_haplotypes():
    _generate_fake_haps()
    _generate_fake_vars()

    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_test_data.hap", pgen=DATADIR / "valhap_test_data.pvar"
        )
        == True
    )


def test_with_empty_lines():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_empty_lines.hap",
        )
        == True
    )


def test_with_out_of_header_metas_sorted():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_out_of_header_metas.hap", sorted=True
        )
        == False
    )


def test_with_out_of_header_metas_unsorted():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_out_of_header_metas.hap", sorted=False
        )
        == True
    )


def test_with_10_extras_reordered():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_10_extras_reordered.hap"
        )
        == True
    )


def test_with_unexistent_reorders():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_unexistent_reorders.hap"
        )
        == False
    )


def test_with_unexistent_fields():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_unexistent_fields.hap"
        )
        == False
    )


def test_with_inadequate_version():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_inadequate_version.hap"
        )
        == False
    )


def test_with_no_version():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_no_version.hap"
        )
        == False
    )


def test_unreadable_hapfile():
    assert val_hapfile.is_hapfile_valid(Path("NON_EXISTENT_FILENAME.hap")) == False
