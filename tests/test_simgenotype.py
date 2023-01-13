from pathlib import Path

import pytest
from click.testing import CliRunner

from haptools.data import Data
from haptools.__main__ import main


# @pytest.mark.skip(reason="this test takes a long time (~2 mins) to run")
def test_basic(capfd):
    prefix = Path("tests/data/example_simgenotype.vcf")

    cmd = " ".join(
        [
            "simgenotype",
            "--model tests/data/outvcf_gen.dat",
            "--mapdir tests/data/map/",
            "--region 1:1-83000",
            "--ref_vcf tests/data/outvcf_test.vcf.gz",
            "--sample_info tests/data/outvcf_info.tab",
            "--pop_field",
            f"--out {prefix}",
        ]
    )
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert result.exit_code == 0
    assert prefix.with_suffix(".bp").exists()
    assert prefix.exists()

    # delete the files and directory we just created
    prefix.with_suffix(".bp").unlink()
    prefix.unlink()


def test_pgen_output(capfd):
    pytest.importorskip("pgenlib")
    prefix = Path("tests/data/example_simgenotype.pgen")

    cmd = " ".join(
        [
            "simgenotype",
            "--model tests/data/outvcf_gen.dat",
            "--mapdir tests/data/map/",
            "--region 1:1-83000",
            "--ref_vcf tests/data/outvcf_test.pgen",
            "--sample_info tests/data/outvcf_info.tab",
            f"--out {prefix}",
        ]
    )
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert result.exit_code == 0
    assert prefix.with_suffix(".bp").exists()
    assert prefix.exists()
    assert prefix.with_suffix(".pvar").exists()
    assert prefix.with_suffix(".psam").exists()

    # delete the files and directory we just created
    prefix.with_suffix(".bp").unlink()
    prefix.unlink()
    prefix.with_suffix(".pvar").unlink()
    prefix.with_suffix(".psam").unlink()
