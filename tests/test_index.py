import shutil
from pathlib import Path
from fileinput import hook_compressed

import numpy as np
from click.testing import CliRunner

from haptools.__main__ import main

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic(capfd):
    file = Path("tests/data/basic.hap")
    tmp_file = Path("test.hap")

    # copy the file so that we don't affect anything in the tests/data directory
    shutil.copy(str(file), str(tmp_file))

    cmd = f"index {tmp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .hap.gz file is the same as the file in tests/data/
    with hook_compressed(tmp_file.with_suffix(".hap.gz"), mode="rt") as haps:
        haps = filter(lambda l: not l.startswith("#"), haps.read().splitlines())
        with hook_compressed(file.with_suffix(".hap.gz"), mode="rt") as expected:
            exp = filter(lambda l: not l.startswith("#"), expected.read().splitlines())
            assert list(haps) == list(exp)
    # check that the output .hap.gz.tbi file is the same, too
    tbi_nocomment = ".nocomment.hap.gz.tbi"
    with hook_compressed(tmp_file.with_suffix(".hap.gz.tbi"), mode="rb") as haps:
        with hook_compressed(file.with_suffix(tbi_nocomment), mode="rb") as expected:
            assert haps.read() == expected.read()
    assert result.exit_code == 0

    tmp_file.unlink()
    tmp_file.with_suffix(".hap.gz").unlink()
    tmp_file.with_suffix(".hap.gz.tbi").unlink()


def test_no_sort(capfd):
    cmd = f"index --no-sort --output test.hap.gz tests/data/basic.hap.gz"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .hap.gz file is the same as the file in tests/data/
    with hook_compressed("test.hap.gz", mode="rt") as haps:
        with hook_compressed("tests/data/basic.hap.gz", mode="rt") as expected:
            assert haps.read() == expected.read()
    # check that the output .hap.gz.tbi file is the same, too
    tbi_nocomment = ".comment.hap.gz.tbi"
    with hook_compressed("test.hap.gz.tbi", mode="rb") as haps:
        with hook_compressed("tests/data/basic" + tbi_nocomment, mode="rb") as expected:
            assert haps.read() == expected.read()
    assert result.exit_code == 0

    Path("test.hap.gz").unlink()
    Path("test.hap.gz").with_suffix(".gz.tbi").unlink()
