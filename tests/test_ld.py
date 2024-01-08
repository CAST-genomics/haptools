from pathlib import Path

from click.testing import CliRunner

from haptools.data import Data
from haptools.__main__ import main

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic(capfd):
    expected = """#\torderH\tld
#\tversion\t0.2.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t21\t26938353\t26938989\tchr21.q.3365*11\t0.995
H\t21\t26938989\t26941960\tchr21.q.3365*10\t-0.012
V\tchr21.q.3365*10\t26938989\t26938989\t21_26938989_G_A\tA
V\tchr21.q.3365*10\t26940815\t26940815\t21_26940815_T_C\tT
V\tchr21.q.3365*10\t26941960\t26941960\t21_26941960_A_G\tA
V\tchr21.q.3365*11\t26938353\t26938353\t21_26938353_T_C\tT
V\tchr21.q.3365*11\t26938989\t26938989\t21_26938989_G_A\tA
"""
    gt_file = DATADIR / "example.vcf.gz"
    hp_file = DATADIR / "basic.hap.gz"

    cmd = f"ld chr21.q.3365*1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0


def test_simple_with_repeat(capfd):
    gt_file = DATADIR / "simple.vcf"
    hp_file = DATADIR / "simple.hap"

    cmd = f"ld H1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out
    expected = captured.out
    assert result.exit_code == 0

    # now try again. The result should be the same because it ignores repeats
    hp_file = DATADIR / "simple_tr.hap"

    cmd = f"ld H1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert expected == captured.out
    assert result.exit_code == 0


def test_basic_variant(capfd):
    expected = """#\torderH\tld
#\tversion\t0.2.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t19\t45411941\t45412079\tAPOe4\t0.999
V\tAPOe4\t45411941\t45411941\trs429358\tC
V\tAPOe4\t45412079\t45412079\trs7412\tC
"""
    tmp_file = Path("apoe4_ld.hap")
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld -o {tmp_file} rs429358 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    assert result.exit_code == 0

    with Data.hook_compressed(tmp_file, mode="r") as haps:
        assert haps.read() == expected

    tmp_file.unlink()


def test_from_gts(capfd):
    expected = """CHR\tBP\tSNP\tR
19\t45411941\trs429358\t0.999
19\t45411947\trs11542041\t0.027
19\t45411962\trs573658040\t-0.012
19\t45411965\trs543363163\t-0.012
19\t45412006\trs563140413\t-0.012
19\t45412007\trs531939919\t-0.012
19\t45412040\trs769455\t0.006
19\t45412079\trs7412\t-0.098
"""
    tmp_file = Path("apoe4.ld")
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld --from-gts -o apoe4.ld APOe4 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    assert result.exit_code == 0

    with Data.hook_compressed(tmp_file, mode="r") as snps:
        assert snps.read() == expected

    tmp_file.unlink()


def test_from_gts_ids(capfd):
    expected = """CHR\tBP\tSNP\tR
19\t45411965\trs543363163\t-0.012
19\t45412079\trs7412\t-0.098
"""
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld --from-gts -i rs543363163 -i rs7412 APOe4 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0
