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
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_10_extras_reordered.hap")
        == True
    )


def test_with_unexistent_reorders():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_unexistent_reorders.hap")
        == False
    )


def test_with_unexistent_fields():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_unexistent_fields.hap")
        == False
    )


def test_with_inadequate_version():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_inadequate_version.hap")
        == False
    )


def test_with_no_version():
    assert val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_no_version.hap") == False


def test_with_multiple_versions():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_multiple_versions.hap")
        == False
    )


def test_with_inadequate_version_columns():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_inadequate_version_columns.hap"
        )
        == False
    )


def test_with_invalid_column_addition_column_count():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_invalid_column_addition_column_count.hap"
        )
        == False
    )


def test_with_invalid_column_addition_types():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_invalid_column_addition_types.hap"
        )
        == False
    )


def test_with_invalid_column_addition_data_types():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_invalid_column_addition_data_types.hap"
        )
        == False
    )


def test_with_insufficient_columns():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_insufficient_columns.hap")
        == False
    )


def test_with_inconvertible_starts():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_inconvertible_starts.hap")
        == False
    )


def test_with_inconvertible_ends():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_inconvertible_ends.hap")
        == False
    )


def test_with_inconvertible_starts_var():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_inconvertible_starts_var.hap"
        )
        == False
    )


def test_with_inconvertible_ends_var():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_inconvertible_ends_var.hap")
        == False
    )


def test_valhap_with_start_after_end():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_start_after_end.hap")
        == False
    )


def test_is_directory():
    assert val_hapfile.is_hapfile_valid(DATADIR / "valhap_is_directory.hap") == False


def test_with_variant_id_of_chromosome():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_variant_id_of_chromosome.hap"
        )
        == False
    )


def test_with_hrid_of_chromosome():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_hrid_of_chromosome.hap")
        == False
    )


def test_with_unexistent_col_in_order():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_unexistent_col_in_order.hap"
        )
        == False
    )


def test_with_unassociated_haplotype():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_unassociated_haplotype.hap")
        == False
    )


def test_with_unrecognizable_allele():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_unrecognizable_allele.hap")
        == False
    )


def test_with_duplicate_ids():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_duplicate_ids.hap") == False
    )


def test_with_duplicate_vids_per_haplotype():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_duplicate_vids_per_haplotype.hap"
        )
        == False
    )


def test_with_excol_of_wrong_type():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_excol_of_wrong_type.hap")
        == False
    )


def test_with_multiple_order_defs():
    assert (
        val_hapfile.is_hapfile_valid(DATADIR / "valhap_with_multiple_order_defs.hap")
        == False
    )


def test_with_insufficient_excols_in_reorder():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_insufficient_excols_in_reorder.hap"
        )
        == False
    )


def test_with_variant_inexistent_haplotype_id():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "valhap_with_variant_inexistent_haplotype_id.hap"
        )
        == False
    )


def test_with_missing_variant_in_pvar():
    assert (
        val_hapfile.is_hapfile_valid(
            DATADIR / "simple.hap", pgen=DATADIR / "basic_missing_ids.pvar"
        )
        == False
    )


def test_unreadable_hapfile():
    assert val_hapfile.is_hapfile_valid(Path("NON_EXISTENT_FILENAME.hap")) == False
