import shutil
from pathlib import Path

import numpy as np
from click.testing import CliRunner

from haptools.data import Data
from haptools.__main__ import main

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic(capfd):
    file = DATADIR / "basic.hap"
    tmp_file = Path("test.hap")

    # copy the file so that we don't affect anything in the tests/data directory
    shutil.copy(str(file), str(tmp_file))

    cmd = f"index {tmp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .hap.gz file is the same as the file in tests/data/
    with Data.hook_compressed(tmp_file.with_suffix(".hap.gz"), mode="rt") as haps:
        haps = filter(lambda l: not l.startswith("#"), haps.read().splitlines())
        with Data.hook_compressed(file.with_suffix(".hap.gz"), mode="rt") as expected:
            exp = filter(lambda l: not l.startswith("#"), expected.read().splitlines())
            assert list(haps) == list(exp)
    # check that the output .hap.gz.tbi file exists
    assert tmp_file.with_suffix(".hap.gz.tbi").is_file()
    assert result.exit_code == 0

    tmp_file.unlink()
    tmp_file.with_suffix(".hap.gz").unlink()
    tmp_file.with_suffix(".hap.gz.tbi").unlink()


def test_basic_w_output(capfd):
    file = DATADIR / "basic.hap"
    tmp_file = Path("test.hap")
    tmp_file_out = Path("sorted.test.hap.gz")

    # copy the file so that we don't affect anything in the tests/data directory
    shutil.copy(str(file), str(tmp_file))

    cmd = f"index --output {tmp_file_out} {tmp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .hap.gz file is the same as the file in tests/data/
    with Data.hook_compressed(tmp_file_out, mode="rt") as haps:
        haps = filter(lambda l: not l.startswith("#"), haps.read().splitlines())
        with Data.hook_compressed(file.with_suffix(".hap.gz"), mode="rt") as expected:
            exp = filter(lambda l: not l.startswith("#"), expected.read().splitlines())
            assert list(haps) == list(exp)
    # check that the output .hap.gz.tbi file exists
    assert tmp_file_out.with_suffix(".gz.tbi").is_file()
    assert result.exit_code == 0

    tmp_file.unlink()
    tmp_file_out.unlink()
    tmp_file_out.with_suffix(".gz.tbi").unlink()


def test_no_sort(capfd):
    file = DATADIR / "basic.hap.gz"
    cmd = f"index --no-sort --output test.hap.gz {file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .hap.gz file is the same as the file in tests/data/
    with Data.hook_compressed("test.hap.gz", mode="rt") as haps:
        with Data.hook_compressed(file, mode="rt") as expected:
            assert haps.read() == expected.read()
    # check that the output .hap.gz.tbi exists
    assert Path("test.hap.gz.tbi").is_file()
    assert result.exit_code == 0

    Path("test.hap.gz").unlink()
    Path("test.hap.gz").with_suffix(".gz.tbi").unlink()
