from pathlib import Path

import pytest

from . import test_data
from haptools import validate as val_hapfile

DATADIR = Path(__file__).parent.joinpath("data") / "valhap"


def test_generated_haplotypes():
    datadir = Path(__file__).parent.joinpath("data")
    hapfile = Path(datadir / "simple.hap")
    pvarfile = Path(datadir / "simple.pvar")

    assert val_hapfile.is_hapfile_valid(hapfile, pvar=pvarfile)


def test_with_empty_lines():
    assert val_hapfile.is_hapfile_valid(DATADIR / "empty_lines.hap")


def test_with_out_of_header_metas_sorted():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "out_of_header_metas.hap", sorted=True
    )


def test_with_out_of_header_metas_unsorted():
    assert val_hapfile.is_hapfile_valid(
        DATADIR / "out_of_header_metas.hap", sorted=False
    )


def test_with_10_extras_reordered():
    assert val_hapfile.is_hapfile_valid(DATADIR / "10_extras_reordered.hap")


def test_with_unexistent_reorders():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "unexistent_reorders.hap")


def test_with_unexistent_fields():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "unexistent_fields.hap")


def test_with_inadequate_version():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inadequate_version.hap")


def test_with_no_version():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "no_version.hap")


def test_with_multiple_versions():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "multiple_versions.hap")


def test_with_inadequate_version_columns():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inadequate_version_columns.hap")


def test_with_invalid_column_addition_column_count():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "invalid_column_addition_column_count.hap"
    )


def test_with_invalid_column_addition_types():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "invalid_column_addition_types.hap"
    )


def test_with_invalid_column_addition_data_types():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "invalid_column_addition_data_types.hap"
    )


def test_with_insufficient_columns():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "insufficient_columns.hap")


def test_with_inconvertible_starts():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inconvertible_starts.hap")


def test_with_inconvertible_ends():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inconvertible_ends.hap")


def test_with_inconvertible_starts_var():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inconvertible_starts_var.hap")


def test_with_inconvertible_ends_var():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "inconvertible_ends_var.hap")


def test_start_after_end():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "start_after_end.hap")


def test_is_directory():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "is_directory.hap")


def test_with_variant_id_of_chromosome():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "variant_id_of_chromosome.hap")


def test_with_hrid_of_chromosome():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "hrid_of_chromosome.hap")


def test_with_unexistent_col_in_order():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "unexistent_col_in_order.hap")


def test_with_unassociated_haplotype():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "unassociated_haplotype.hap")


def test_with_unrecognizable_allele():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "unrecognizable_allele.hap")


def test_with_duplicate_ids():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "duplicate_ids.hap")


def test_with_duplicate_vids_per_haplotype():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "duplicate_vids_per_haplotype.hap"
    )


def test_with_excol_of_wrong_type():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "excol_of_wrong_type.hap")


def test_with_multiple_order_defs():
    assert not val_hapfile.is_hapfile_valid(DATADIR / "multiple_order_defs.hap")


def test_with_insufficient_excols_in_reorder():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "insufficient_excols_in_reorder.hap"
    )


def test_with_variant_inexistent_haplotype_id():
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "variant_inexistent_haplotype_id.hap"
    )


def test_with_missing_variant_in_pvar():
    pgenlib = pytest.importorskip("pgenlib")
    assert not val_hapfile.is_hapfile_valid(
        DATADIR / "simple.hap", pvar=DATADIR / "basic_missing_ids.pvar"
    )


def test_unreadable_hapfile():
    assert not val_hapfile.is_hapfile_valid(Path("NON_EXISTENT_FILENAME.hap"))
