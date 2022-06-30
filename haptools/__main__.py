#!/usr/bin/env python

from __future__ import annotations
import sys
import click
from pathlib import Path

# AVOID IMPORTING ANYTHING HERE
# any imports we put here will make it slower to use the command line client
# a basic "haptools --help" should be quick and require very few imports, for example

################### Haptools ##################
@click.group()
@click.version_option()
def main():
    """
    haptools: A toolkit for simulating and analyzing genotypes and 
    phenotypes while taking into account haplotype information
    """
    pass

############ Haptools karyogram ###############
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
    from .karyogram import PlotKaryogram
    if colors is not None:
        colors = dict([item.split(":") for item in colors.split(",")])
    PlotKaryogram(bp, sample, out, \
        centromeres_file=centromeres, title=title, colors=colors)

############ Haptools simgenotype ###############
@main.command()
@click.option('--model', help="Admixture model in .dat format. See docs for info.", \
    type=str, required=True)
@click.option('--mapdir', help="Directory containing files with chr\{1-22,X\} and ending in .map in the file name with genetic map coords.", \
    required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True))
@click.option('--out', help="Prefix to name output files.", \
    type=str, required=True)
@click.option('--chroms', help='Sorted and comma delimited list of chromosomes to simulate. ex: 1,2,3,5,6,21,X', \
    type=str, default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X", required=True)
@click.option('--seed', help="Random seed. Set to make simulations reproducible", \
    type=int, required=False, default=None)
@click.option('--popsize', help="Number of samples to simulate each generation", \
    type=int, required=False, default=10000, hidden=True)
@click.option('--invcf', help="VCF file used as reference for creation of simulated samples respective genotypes.", required=True)
@click.option('--sample_info', help="File that maps samples from the reference VCF (--invcf) to population codes " \
              "describing the populations in the header of the model file.", required=True)
def simgenotype(invcf, sample_info, model, mapdir, out, popsize, seed, chroms):
    """
    Simulate admixed genomes under a pre-defined model.

    Example:
    haptools simgenotype \
      --model ./tests/data/outvcf_gen.dat \
      --mapdir ./tests/map/ \
      --chroms 1,2 \
      --invcf ./tests/data/outvcf_test.vcf \
      --sample_info ./tests/data/outvcf_info.tab \
      --out ./tests/data/example_simgenotype
    """
    from .sim_genotype import simulate_gt, write_breakpoints, output_vcf, validate_params
    chroms = chroms.split(',')
    validate_params(model, mapdir, chroms, popsize, invcf, sample_info)
    samples, breakpoints = simulate_gt(model, mapdir, chroms, popsize, seed)
    breakpoints = write_breakpoints(samples, breakpoints, out)
    output_vcf(breakpoints, model, invcf, sample_info, out)

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
    from .sim_phenotypes import simulate_pt
    # Basic checks on input
    # TODO - check VCF zipped, check only one of simu-qt/simu-cc,
    # check values of other inputs
    # Only use simu-k for case/control

    # Run simulation
    simulate_pt(vcf, hap, simu_rep, \
        simu_hsq, simu_k, simu_qt, simu_cc, out)

@main.command(short_help="Transform a genotypes matrix via a set of haplotypes")
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help="""
    The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'\n
    For this to work, the VCF must be indexed and the seqname must match!""",
)
@click.option(
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A VCF file containing haplotype 'genotypes'",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default="only errors",
    help="The level of verbosity desired",
)
def transform(
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    output: Path = Path("-"),
    verbosity: str = 'CRITICAL',
):
    """
    Creates a VCF composed of haplotypes

    GENOTYPES must be formatted as a VCF and HAPLOTYPES must be formatted according
    to the .hap format spec

    \f
    Examples
    --------
    >>> haptools transform tests/data/example.vcf.gz tests/data/example.hap.gz > example_haps.vcf

    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF format
    haplotypes : Path
        The path to the haplotypes in a .hap file
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    sample : Tuple[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    samples_file : Path, optional
        A single column txt file containing a list of the samples (one per line) to
        subset from the genotypes file
    output : Path, optional
        The location to which to write output
    verbosity : str, optional
        The level of verbosity desired in messages written to stderr
    """
    import logging

    from haptools import data
    from .haplotype import HaptoolsHaplotype

    log = logging.getLogger("run")
    logging.basicConfig(
        format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
    )
    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = samps_file.read().splitlines()
    elif samples:
        # needs to be converted from tuple to list
        samples = list(samples)
    else:
        samples = None

    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes)
    hp.read(region=region)
    log.info("Extracting variants from haplotypes")
    variants = {var.id for hap in hp.data.values() for var in hap.variants}
    log.info("Loading genotypes")
    gt = data.GenotypesRefAlt(genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing(discard_also=True)
    gt.check_biallelic(discard_also=True)
    gt.check_phase()
    log.info("Transforming genotypes via haplotypes")
    hp_gt = data.GenotypesRefAlt(fname=output, log=log)
    hp.transform(gt, hp_gt)
    log.info("Writing haplotypes to VCF")
    hp_gt.write()


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")
