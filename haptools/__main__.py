#!/usr/bin/env python

import sys
import click

from pathlib import Path
from typing import Union, Tuple

#from . import data, tree
sys.path.append('./admixture_sim')
from .admixture_sim.sim_admixture import simulate_gt, write_breakpoint


@click.group()
@click.version_option()
def main():
    """
    simwas: Simulate phenotypes for GWAS and subsequent fine-mapping

    Use real variants to simulate real, biological LD patterns.
    """
    pass


@main.command()
@click.option('--invcf')
@click.option('--sample_info')
@click.option('--model', required=True)
@click.option('--coords', required=True)
@click.option('--out', required=True)
def simulate():
    """
    Use the tool to simulate genotypes
    """
    breakpoints = simulate_gt(model_file, coords_file, seed=None)
    write_breakpoint(breakpoints, out)


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m simwas' on the command line
    main(prog_name="simwas")
