from pathlib import Path

import pytest
from click.testing import CliRunner

from haptools.__main__ import main
from haptools.validate import is_hapfile_valid

PARENT_DATADIR = Path(__file__).parent.joinpath("data")
DATADIR = Path(__file__).parent.joinpath("data") / "valhap"


def test_generated_haplotypes():
    hapfile = Path(PARENT_DATADIR / "simple.hap")
    pvarfile = Path(PARENT_DATADIR / "simple.pvar")

    assert is_hapfile_valid(hapfile, pvar=pvarfile)


def test_with_empty_lines():
    assert is_hapfile_valid(DATADIR / "empty_lines.hap")


def test_with_out_of_header_metas_sorted():
    assert not is_hapfile_valid(DATADIR / "out_of_header_metas.hap", sorted=True)


def test_with_out_of_header_metas_unsorted():
    assert is_hapfile_valid(DATADIR / "out_of_header_metas.hap", sorted=False)


def test_with_10_extras_reordered():
    assert is_hapfile_valid(DATADIR / "10_extras_reordered.hap")


def test_with_unexistent_reorders():
    assert not is_hapfile_valid(DATADIR / "unexistent_reorders.hap")


def test_with_unexistent_fields():
    assert not is_hapfile_valid(DATADIR / "unexistent_fields.hap")


def test_with_inadequate_version():
    assert not is_hapfile_valid(DATADIR / "inadequate_version.hap")


def test_with_no_version():
    assert not is_hapfile_valid(DATADIR / "no_version.hap")


def test_with_multiple_versions():
    assert not is_hapfile_valid(DATADIR / "multiple_versions.hap")


def test_with_inadequate_version_columns():
    assert not is_hapfile_valid(DATADIR / "inadequate_version_columns.hap")


def test_with_invalid_column_addition_column_count():
    assert not is_hapfile_valid(DATADIR / "invalid_column_addition_column_count.hap")


def test_with_invalid_column_addition_types():
    assert not is_hapfile_valid(DATADIR / "invalid_column_addition_types.hap")


def test_with_invalid_column_addition_data_types():
    assert not is_hapfile_valid(DATADIR / "invalid_column_addition_data_types.hap")


def test_with_insufficient_columns():
    assert not is_hapfile_valid(DATADIR / "insufficient_columns.hap")


def test_with_inconvertible_starts():
    assert not is_hapfile_valid(DATADIR / "inconvertible_starts.hap")


def test_with_inconvertible_ends():
    assert not is_hapfile_valid(DATADIR / "inconvertible_ends.hap")


def test_with_inconvertible_starts_var():
    assert not is_hapfile_valid(DATADIR / "inconvertible_starts_var.hap")


def test_with_inconvertible_ends_var():
    assert not is_hapfile_valid(DATADIR / "inconvertible_ends_var.hap")


def test_start_after_end():
    assert not is_hapfile_valid(DATADIR / "start_after_end.hap")


def test_is_directory():
    assert not is_hapfile_valid(DATADIR / "is_directory.hap")


def test_with_variant_id_of_chromosome():
    assert not is_hapfile_valid(DATADIR / "variant_id_of_chromosome.hap")


def test_with_hrid_of_chromosome():
    assert not is_hapfile_valid(DATADIR / "hrid_of_chromosome.hap")


def test_with_unexistent_col_in_order():
    assert not is_hapfile_valid(DATADIR / "unexistent_col_in_order.hap")


def test_with_unassociated_haplotype():
    assert not is_hapfile_valid(DATADIR / "unassociated_haplotype.hap")


def test_with_unrecognizable_allele():
    assert not is_hapfile_valid(DATADIR / "unrecognizable_allele.hap")


def test_with_duplicate_ids():
    assert not is_hapfile_valid(DATADIR / "duplicate_ids.hap")


def test_with_duplicate_vids_per_haplotype():
    assert not is_hapfile_valid(DATADIR / "duplicate_vids_per_haplotype.hap")


def test_with_excol_of_wrong_type():
    assert not is_hapfile_valid(DATADIR / "excol_of_wrong_type.hap")


def test_with_multiple_order_defs():
    assert not is_hapfile_valid(DATADIR / "multiple_order_defs.hap")


def test_with_insufficient_excols_in_reorder():
    assert not is_hapfile_valid(DATADIR / "insufficient_excols_in_reorder.hap")


def test_with_variant_inexistent_haplotype_id():
    assert not is_hapfile_valid(DATADIR / "variant_inexistent_haplotype_id.hap")


def test_with_missing_variant_in_pvar():
    pgenlib = pytest.importorskip("pgenlib")
    assert not is_hapfile_valid(
        DATADIR / "simple.hap", pvar=DATADIR / "basic_missing_ids.pvar"
    )


def test_unreadable_hapfile():
    assert not is_hapfile_valid(Path("NON_EXISTENT_FILENAME.hap"))


def test_basic(capfd):
    hp_file = DATADIR / "basic.hap"

    cmd = f"validate {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code == 0


def test_no_version(capfd):
    hp_file = DATADIR / "no_version.hap"

    cmd = f"validate {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code != 0


def test_no_version(capfd):
    hp_file = DATADIR / "no_version.hap"

    cmd = f"validate {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code != 0


def test_sorted(capfd):
    hp_file = DATADIR / "out_of_header_metas.hap"

    cmd = f"validate --sorted {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code != 0

    cmd = f"validate --no-sorted {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code == 0


def test_with_pvar(capfd):
    gt_file = PARENT_DATADIR / "simple.pvar"
    hp_file = PARENT_DATADIR / "simple.hap"

    cmd = f"validate --genotypes {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    assert result.exit_code == 0
