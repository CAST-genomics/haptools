#!/usr/bin/env python

"""
TODO:
- add defaults to help messages
- clean up help messages
- group options
"""

import sys
import click

from .sim_admixture import simulate_gt, write_breakpoints
from .karyogram import PlotKaryogram
from .sim_phenotypes import simulate_pt

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
@click.option('--popsize', default=10000, hidden=True)
@click.option('--seed', default=None)
@click.option('--chroms', help='Sorted (1-22, X) and comma delimited list of chromosomes used to simulate admixture. ex: 1,2,3,5,6,21,X \
                                .', type=str, required=True)
def simgenotype(invcf, sample_info, model, mapdir, out, popsize, seed, chroms):
    """
    Use the tool to simulate genotypes
    """
    chroms = chroms.split(',')
    samples, breakpoints = simulate_gt(model, mapdir, chroms, popsize, seed)
    breakpoints = write_breakpoints(samples, breakpoints, out)
    

@main.command()
@click.option('--bp', help="Path to .bp file with breakpoints", \
    type=str, required=True)
@click.option('--sample', help="Sample ID to plot", \
    type=str, required=True)
@click.option('--out', help="Name of output file", \
    type=str, required=True)
@click.option('--title', help="Optional plot title", \
    type=str, required=False)
@click.option('--centromeres', help="Optional file with telomere/centromere cM positions", \
    type=str, required=False)
@click.option('--colors', help="Optional color dictionary. Format is e.g. 'YRI:blue,CEU:green'", \
    type=str, required=False)
def karyogram(bp, sample, out, title, centromeres, colors):
    """
    Visualize a karyogram of local ancestry tracks

    Example:
    haptools karyogram --bp tests/data/5gen.bp --sample Sample_1 \
       --out test.png --centromeres tests/data/centromeres_hg19.txt \
       --colors 'CEU:blue,YRI:red'
    """
    if colors is not None:
        colors = dict([item.split(":") for item in colors.split(",")])
    PlotKaryogram(bp, sample, out, \
        centromeres_file=centromeres, title=title, colors=colors)

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
