#!/usr/bin/env python

import sys
import click

from pathlib import Path
#from typing import Union, Tuple

#from . import data, tree
sys.path.append('./admixture_sim')
from .simulate.sim_admixture import simulate_gt, write_breakpoints


@click.group()
@click.version_option()
def main():
    """
    haptools: Simulate phenotypes for GWAS and subsequent fine-mapping

    Use real variants to simulate real, biological LD patterns.
    """
    pass

@main.command()
@click.option('--invcf')
@click.option('--sample_info')
@click.option('--model', required=True)
@click.option('--coords', required=True)
@click.option('--out', required=True)
def simulate(invcf, sample_info, model, coords, out):
    """
    Use the tool to simulate genotypes
    """
    breakpoints = simulate_gt(model, coords)
    write_breakpoint(breakpoints, out)


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")




