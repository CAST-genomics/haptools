#!/usr/bin/env python

import sys
import click

from pathlib import Path
#from typing import Union, Tuple

#from . import data, tree
sys.path.append('./admixture_sim')
from .simgenotype.sim_admixture import simulate_gt, write_breakpoints
from .karyogram.karyogram import plot_karyogram


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
@click.option('--mapdir', required=True)
@click.option('--out', required=True)
def simgenotype(invcf, sample_info, model, mapdir, out):
    """
    Use the tool to simulate genotypes
    """
    samples, breakpoints = simulate_gt(model, mapdir)
    write_breakpoints(samples, breakpoints, out)

@main.command()
@click.option('--sample_name', type=str)
@click.option('--chrX', default=False, type=bool)
@click.option('--sample_file', required=True)
@click.option('--title', required=True)
@click.option('--centromeres', required=True)
@click.option('--out', required=True)
def karyogram(sample_name, chrx, sample_file, title, centromeres, out):
    """
    Use the tool to visualize breakpoints.
    """
    plot_karyogram(sample_file, title, centromeres, out)
    #plot_karyogram(sample_file, title, centromeres, out, sample_name="Sample_1" chrX=False, colors=None)

if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")




