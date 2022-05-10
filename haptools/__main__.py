#!/usr/bin/env python

"""
TODO:
- add defaults to help messages
- clean up help messages
- group options
"""

from __future__ import annotations
import sys
import click
from pathlib import Path

# AVOID IMPORTING ANYTHING HERE
# any imports you put here will make it slower to use the command line client
# a basic haptools --help should be quick and require very few imports, for example


@click.group()
@click.version_option()
def main():
    """
    haptools: Simulate genotypes and phenotypes for GWAS and subsequent fine-mapping

    Use real variants to simulate biological LD patterns and traits.
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
    from .sim_admixture import simulate_gt, write_breakpoints
    samples, breakpoints = simulate_gt(model, mapdir)
    write_breakpoints(samples, breakpoints, out)

@main.command()
@click.option('--sample_name', type=str)
@click.option('--chrX', is_flag=True)
@click.option('--sample_file', required=True)
@click.option('--title', required=True)
@click.option('--centromeres', required=True)
@click.option('--out', required=True)
def karyogram(sample_name, chrx, sample_file, title, centromeres, out):
    """
    Use the tool to visualize breakpoints.
    """
    from .karyogram import plot_karyogram
    plot_karyogram(sample_file, title, centromeres, out, sample_name, chrx)

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

@main.command()
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
    type=click.File("w"),
    default="-",
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
    output: TextIO = sys.stdout,
    verbosity: str = 'CRITICAL',
):
    """
    Transform a genotypes matrix via a set of haplotypes

    GENOTYPES must be formatted as a VCF and

    HAPLOTYPES must be formatted according to the .hap format spec

    Ex: haptools transform tests/data/example.vcf.gz tests/data/example.hap.gz > example_haps.vcf
    \f
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
    output : TextIO, optional
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
    # load data
    log.info("Loading genotypes")
    gt = data.Genotypes(genotypes)
    gt.read(region=region, samples=samples)
    log.info("Discarding multiallelic variants")
    gt.check_biallelic(discard_also=True)
    gt.check_phase()
    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes, haplotype=HaptoolsHaplotype)
    hp.read(region=region)
    hp_gt = hp.transform(gt)
    # TODO: write hp_gt to output


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")
