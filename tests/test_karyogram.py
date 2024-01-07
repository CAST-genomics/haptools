from pathlib import Path

from click.testing import CliRunner

from haptools.__main__ import main
from haptools.karyogram import GetHaplotypeBlocks, PlotKaryogram

DATADIR = Path(__file__).parent.joinpath("data")


def test_GetHaplotypeBlocks():
    test_file = DATADIR.joinpath("test.bp")
    sample_blocks = GetHaplotypeBlocks(test_file, "Sample_1")
    assert sample_blocks[0][0]["pop"] == "YRI"
    assert sample_blocks[0][0]["chrom"] == 1
    assert sample_blocks[0][0]["start"] == 0.0001
    assert sample_blocks[0][0]["end"] == 168.003442

    assert sample_blocks[1][0]["pop"] == "YRI"
    assert sample_blocks[1][0]["chrom"] == 1
    assert sample_blocks[1][0]["start"] == 0.0001
    assert sample_blocks[1][0]["end"] == 87.107755

    sample_blocks = GetHaplotypeBlocks(test_file, "Sample_2")
    assert sample_blocks[0][-1]["pop"] == "YRI"
    assert sample_blocks[0][-1]["chrom"] == 2
    assert sample_blocks[0][-1]["start"] == 180.837755 + 0.0001
    assert sample_blocks[0][-1]["end"] == 244.341689

    assert sample_blocks[1][0]["pop"] == "YRI"
    assert sample_blocks[1][0]["chrom"] == 1
    assert sample_blocks[1][0]["start"] == 0.0001
    assert sample_blocks[1][0]["end"] == 85.107755


def test_centromere_chroms(capfd):
    tmp_file = Path("test_karyogram.png")

    cmd = " ".join(
        [
            "karyogram",
            "--bp tests/data/outvcf_test.bp",
            "--sample Sample_1",
            f"--out {tmp_file}",
            "--centromeres tests/data/centromeres_hg19.txt",
            "--title Non_23_Chrom_Karyogram",
            "--colors CEU:blue,YRI:red",
        ]
    )
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert result.exit_code == 0

    # delete the file we just created
    tmp_file.unlink()


def test_basic(capfd):
    tmp_file = Path("test_karyogram.png")

    cmd = " ".join(
        [
            "karyogram",
            "--bp tests/data/5gen.bp",
            "--sample Sample_1",
            f"--out {tmp_file}",
            "--centromeres tests/data/centromeres_hg19.txt",
            "--title 5_Generation_Karyogram",
            "--colors CEU:blue,YRI:red",
        ]
    )
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert result.exit_code == 0

    # delete the file we just created
    tmp_file.unlink()
