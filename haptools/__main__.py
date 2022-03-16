#!/usr/bin/env python

"""
TODO:
- add defaults to help messages
- clean up help messages
- group options
- add tabix to haptools-dev?
"""

import sys
import click

from .simgenotype.sim_admixture import simulate_gt, write_breakpoints
from .karyogram.karyogram import plot_karyogram
from .simphenotype.sim_phenotypes import simulate_pt

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

############ Haptools simphenotype ###############
DEFAULT_SIMU_REP = 1
DEFAULT_SIMU_HSQ = 0.1
DEFAULT_SIMU_K = 0.1
##################################################
@main.command()
@click.option('--vcf', help='Phased VCF file', type=str, required=True)
@click.option('--hap', help='Haplotype file with effect sizes', \
        type=str, required=True)
@click.option('--out', help='Prefix for output files', \
        type=str, required=True)
@click.option('--simu-qt', help='Simulate a quantitative trait', \
        default=False, is_flag=True)
@click.option('--simu-cc', help='Simulate a case/control trait', \
        default=False, is_flag=True)
@click.option('--simu-rep', help='Number of rounds of simulation to perform', \
        type=int, default=DEFAULT_SIMU_REP)
@click.option('--simu-hsq', help='Trait heritability', \
        type=float, default=DEFAULT_SIMU_HSQ)
@click.option('--simu-k', help='Specify the disease prevalence', \
        type=float, default=DEFAULT_SIMU_K)
def simphenotype(vcf, hap, simu_rep, simu_hsq, simu_k, simu_qt, simu_cc, out):
    """
    Haplotype-aware phenotype simulation
    """
    # Basic checks on input
    # TODO - check VCF zipped, check only one of simu-qt/simu-cc,
    # check values of other inputs
    # Only use simu-k for case/control

    # Run simulation
    simulate_pt(vcf, hap, simu_rep, \
        simu_hsq, simu_k, simu_qt, simu_cc, out)

if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")