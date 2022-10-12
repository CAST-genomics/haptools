from pathlib import Path

import numpy as np
from click.testing import CliRunner

from haptools.__main__ import main

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic(capfd):
    cmd = "ld chr21.q.3365*1 tests/data/example.vcf.gz tests/data/basic.hap.gz"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == """#\torderH\tld
#\tversion\t0.1.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t21\t26938353\t26938989\tchr21.q.3365*11\t0.995
H\t21\t26938989\t26941960\tchr21.q.3365*10\t-0.012
V\tchr21.q.3365*10\t26938989\t26938989\t21_26938989_G_A\tA
V\tchr21.q.3365*10\t26940815\t26940815\t21_26940815_T_C\tT
V\tchr21.q.3365*10\t26941960\t26941960\t21_26941960_A_G\tA
V\tchr21.q.3365*11\t26938353\t26938353\t21_26938353_T_C\tT
V\tchr21.q.3365*11\t26938989\t26938989\t21_26938989_G_A\tA
"""
    assert result.exit_code == 0

def test_basic_variant(capfd):
    cmd = "ld rs429358 tests/data/apoe.vcf.gz tests/data/apoe4.hap"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == """#\torderH\tld
#\tversion\t0.1.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t19\t45411941\t45412079\tAPOe4\t0.999
V\tAPOe4\t45411941\t45411941\trs429358\tC
V\tAPOe4\t45412079\t45412079\trs7412\tC
"""
    assert result.exit_code == 0

def test_from_gts(capfd):
    cmd = "ld --from-gts APOe4 tests/data/apoe.vcf.gz tests/data/apoe4.hap"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == """CHR\tBP\tSNP\tR
19\t45411941\trs429358\t0.999
19\t45411947\trs11542041\t0.027
19\t45411962\trs573658040\t-0.012
19\t45411965\trs543363163\t-0.012
19\t45412006\trs563140413\t-0.012
19\t45412007\trs531939919\t-0.012
19\t45412040\trs769455\t0.006
19\t45412079\trs7412\t-0.098
"""
    assert result.exit_code == 0

def test_from_gts_ids(capfd):
    cmd = "ld --from-gts -i rs543363163 -i rs7412 APOe4 tests/data/apoe.vcf.gz tests/data/apoe4.hap"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert captured.out == """CHR\tBP\tSNP\tR
19\t45411965\trs543363163\t-0.012
19\t45412079\trs7412\t-0.098
"""
    assert result.exit_code == 0
