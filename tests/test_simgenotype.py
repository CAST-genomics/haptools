from pathlib import Path

import pytest
from click.testing import CliRunner

from haptools.data import Data
from haptools.__main__ import main

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic(capfd):
    prefix = DATADIR / "example_simgenotype.vcf"
    dat_file = DATADIR / "outvcf_gen.dat"
    map_dir = DATADIR / "map"
    ref_vcf_file = DATADIR / "outvcf_test.vcf.gz"
    samp_info_file = DATADIR / "outvcf_info.tab"

    cmd = " ".join(
        [
            "simgenotype",
            f"--model {dat_file}",
            f"--mapdir {map_dir}",
            "--region 1:1-83000",
            f"--ref_vcf {ref_vcf_file}",
            f"--sample_info {samp_info_file}",
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
    prefix = DATADIR / "example_simgenotype.pgen"
    dat_file = DATADIR / "outvcf_gen.dat"
    map_dir = DATADIR / "map"
    ref_vcf_file = DATADIR / "outvcf_test.pgen"
    samp_info_file = DATADIR / "outvcf_info.tab"

    cmd = " ".join(
        [
            "simgenotype",
            f"--model {dat_file}",
            f"--mapdir {map_dir}",
            "--region 1:1-83000",
            f"--ref_vcf {ref_vcf_file}",
            f"--sample_info {samp_info_file}",
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


def test_pgen_output_var_greater(capfd):
    prefix = DATADIR / "example_simgenotype.pgen"
    dat_file = DATADIR / "outvcf_gen.dat"
    map_dir = DATADIR / "map"
    ref_vcf_file = DATADIR / "var_greater.vcf.gz"
    samp_info_file = DATADIR / "outvcf_info.tab"

    cmd = " ".join(
        [
            "simgenotype",
            f"--model {dat_file}",
            f"--mapdir {map_dir}",
            "--region 1:249320800-249403800",
            f"--ref_vcf {ref_vcf_file}",
            f"--sample_info {samp_info_file}",
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


def test_pgen_output_chunked(capfd):
    prefix = DATADIR / "example_simgenotype.pgen"
    dat_file = DATADIR / "outvcf_gen.dat"
    map_dir = DATADIR / "map"
    ref_vcf_file = DATADIR / "outvcf_test.pgen"
    samp_info_file = DATADIR / "outvcf_info.tab"

    cmd = " ".join(
        [
            "simgenotype",
            f"--model {dat_file}",
            f"--mapdir {map_dir}",
            "--region 1:1-83000",
            f"--ref_vcf {ref_vcf_file}",
            f"--sample_info {samp_info_file}",
            f"--out {prefix}",
            f"--chunk-size 1",
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
