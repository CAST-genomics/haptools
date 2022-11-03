#!/usr/bin/env python

from __future__ import annotations
import sys
from pathlib import Path

import click

# AVOID IMPORTING ANYTHING ABOVE
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


@main.command()
@click.option(
    "--bp",
    type=str,
    required=True,
    help="Path to .bp file with breakpoints",
)
@click.option(
    "--sample",
    type=str,
    required=True,
    help="Sample ID to plot",
)
@click.option(
    "--out",
    type=str,
    required=True,
    help="Name of output file",
)
@click.option(
    "--title",
    type=str,
    required=False,
    help="Optional plot title",
)
@click.option(
    "--centromeres",
    type=str,
    required=False,
    help="Optional file with telomere/centromere cM positions",
)
@click.option(
    "--colors",
    type=str,
    required=False,
    help="Optional color dictionary. Format is e.g. 'YRI:blue,CEU:green'",
)
def karyogram(bp, sample, out, title, centromeres, colors):
    """
    Visualize a karyogram of local ancestry tracks
    """
    from .karyogram import PlotKaryogram

    if colors is not None:
        colors = dict([item.split(":") for item in colors.split(",")])
    PlotKaryogram(
        bp, sample, out, centromeres_file=centromeres, title=title, colors=colors
    )


@main.command()
@click.option(
    "--model",
    type=str,
    required=True,
    help="Admixture model in .dat format. See docs for info.",
)
@click.option(
    "--mapdir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help=(
        "Directory containing files with chr{1-22,X} and ending in .map in the file "
        "name with genetic map coords."
    ),
)
@click.option(
    "--out",
    type=str,
    required=True,
    help="Prefix to name output files.",
)
@click.option(
    "--chroms",
    type=str,
    required=False,
    default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
    help="Sorted and comma delimited list of chromosomes to simulate",
)
@click.option(
    "--seed",
    type=int,
    required=False,
    default=None,
    help="Random seed. Set to make simulations reproducible",
)
@click.option(
    "--popsize",
    type=int,
    hidden=True,
    default=10000,
    required=False,
    help="Number of samples to simulate each generation",
)
@click.option(
    "--invcf",
    required=True,
    help=(
        "VCF file used as reference for creation of simulated samples respective "
        "genotypes."
    ),
)
@click.option(
    "--sample_info",
    required=True,
    help=(
        "File that maps samples from the reference VCF (--invcf) to population codes "
        "describing the populations in the header of the model file."
    ),
)
@click.option(
    "--region",
    required=False,
    default=None,
    help=(
        "Subset the simulation to a specific region in a chromosome using the form"
        " chrom:start-end. Example 2:1000-2000"
    ),
)
@click.option(
    "--only_breakpoint",
    hidden=True,
    is_flag=True,
    required=False,
    help=(
        "Flag used to determine whether to only output breakpoints or continue to "
        "simulate a vcf file."
    ),
)
@click.option(
    "--verbose",
    hidden=True,
    is_flag=True,
    required=False,
    help=(
        "Output time metrics for each section, breakpoint simulation, vcf creation, "
        "and total exection."
    ),
)
def simgenotype(
    invcf,
    sample_info,
    model,
    mapdir,
    out,
    popsize,
    seed,
    chroms,
    region,
    only_breakpoint,
    verbose,
):
    """
    Simulate admixed genomes under a pre-defined model.
    """
    import re
    import time
    from .sim_genotype import (
        output_vcf,
        simulate_gt,
        validate_params,
        write_breakpoints,
    )

    start = time.time()

    # parse region and chroms parameters
    if not (chroms or region):
        raise Exception("Either chroms or region must be specified.")
    if region:
        region_info = re.split(":|-", region)
        try:
            region = {
                "chr": region_info[0],
                "start": int(region_info[1]),
                "end": int(region_info[2]),
            }
            chroms = [region["chr"]]
        except:
            raise Exception(
                "Unable to parse region. Please ensure it has the correct format "
                "<chr>:<start>-<end> eg. 1:1-2000"
            )
    else:
        chroms = chroms.split(",")

    # Handle if mapdir has a '/' at the end
    if mapdir[-1] == "/":
        mapdir = mapdir[:-1]

    # simulate breakpoints
    popsize = validate_params(
        model, mapdir, chroms, popsize, invcf, sample_info, region, only_breakpoint
    )
    samples, breakpoints = simulate_gt(model, mapdir, chroms, region, popsize, seed)
    breakpoints = write_breakpoints(samples, breakpoints, out)
    bp_end = time.time()

    # simulate vcfs
    vcf_start = time.time()
    if not only_breakpoint:
        # TODO add region functionality
        output_vcf(breakpoints, chroms, model, invcf, sample_info, region, out)
    end = time.time()

    if verbose:
        print(f"Time elapsed for breakpoint simulation: {bp_end - start}")
        print(f"Time elapse for creating vcf: {end - vcf_start}")
        print(f"Time elapsed for simgenotype execution: {end - start}")


@main.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-r",
    "--replications",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of rounds of simulation to perform",
)
@click.option(
    "-h",
    "--heritability",
    type=click.FloatRange(min=0, max=1),
    default=None,
    show_default=True,
    help="Trait heritability",
)
@click.option(
    "-p",
    "--prevalence",
    type=click.FloatRange(min=0, max=1, min_open=False, max_open=True),
    show_default="quantitative trait",
    help="Disease prevalence if simulating a case-control trait",
)
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
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
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A list of the haplotype IDs from the .hap file to use as causal variables "
        "(ex: '-i H1 -i H2')."
    ),
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype IDs "
        "(one per line) to subset from the .hap file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A TSV file containing simulated phenotypes",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default="only errors",
    help="The level of verbosity desired",
)
def simphenotype(
    genotypes: Path,
    haplotypes: Path,
    replications: int = 1,
    heritability: float = None,
    prevalence: float = None,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    output: Path = Path("-"),
    verbosity: str = "ERROR",
):
    """
    Haplotype-aware phenotype simulation. Create a set of simulated phenotypes from a
    set of haplotypes.

    GENOTYPES must be formatted as a VCF or PGEN file and HAPLOTYPES must be formatted
    according to the .hap format spec

    Note: GENOTYPES must be the output from the transform subcommand.
    """
    import logging

    from .sim_phenotype import simulate_pt

    log = logging.getLogger("haptools simphenotype")
    db_time = "|%(asctime)s" if verbosity == "DEBUG" else ""
    logging.basicConfig(
        format="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
        datefmt="%H:%M:%S",
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

    if ids_file:
        with ids_file as id_file:
            ids = set(id_file.read().splitlines())
    elif ids:
        ids = set(ids)
    else:
        ids = None

    # Run simulation
    simulate_pt(
        genotypes,
        haplotypes,
        replications,
        heritability,
        prevalence,
        region,
        samples,
        ids,
        chunk_size,
        output,
        log,
    )


@main.command(short_help="Transform a genotypes matrix via a set of haplotypes")
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
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
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help="A list of the haplotype IDs to use from the .hap file (ex: '-i H1 -i H2').",
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype IDs "
        "(one per line) to subset from the .hap file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Ignore any samples that are missing genotypes for the required variants",
)
@click.option(
    "--ancestry",
    is_flag=True,
    show_default=True,
    default=False,
    help="Also transform using VCF 'POP' FORMAT field and 'ancestry' .hap extra field",
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
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    ancestry: bool = False,
    output: Path = Path("-"),
    verbosity: str = "CRITICAL",
):
    """
    Creates a VCF composed of haplotypes

    GENOTYPES must be formatted as a VCF or PGEN and HAPLOTYPES must be formatted
    according to the .hap format spec
    """
    import logging

    from .transform import transform_haps

    log = logging.getLogger("haptools transform")
    db_time = "|%(asctime)s" if verbosity == "DEBUG" else ""
    logging.basicConfig(
        format="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
        datefmt="%H:%M:%S",
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

    if ids_file:
        with ids_file as id_file:
            ids = set(id_file.read().splitlines())
    elif ids:
        ids = set(ids)
    else:
        ids = None

    transform_haps(
        genotypes,
        haplotypes,
        region,
        samples,
        ids,
        chunk_size,
        discard_missing,
        ancestry,
        output,
        log,
    )


@main.command(short_help="Compute pair-wise LD")
@click.argument("target", type=str)
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
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
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A list of the haplotype IDs to use from the .hap file (ex: '-i H1 -i H2'). "
        "Or, if --from-gts, a list of the variant IDs to use from the genotypes file."
        "\nFor this to work, the .hap file must be indexed"
    ),
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype (or variant) IDs "
        "(one per line) to subset from the .hap (or genotype) file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Ignore any samples that are missing genotypes for the required variants",
)
@click.option(
    "--from-gts",
    is_flag=True,
    show_default=True,
    default=False,
    help=(
        "By default, LD is computed with the haplotypes in the .hap file. Use this "
        "switch to compute LD with the genotypes in the genotypes file, instead."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .hap file containing haplotypes and their LD with TARGET",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default="only errors",
    help="The level of verbosity desired",
)
def ld(
    target: str,
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    from_gts: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "CRITICAL",
):
    """
    Compute the pair-wise LD (Pearson's correlation) between haplotypes (or variants)
    and a single TARGET haplotype (or variant)

    GENOTYPES must be formatted as a VCF or PGEN and HAPLOTYPES must be formatted
    according to the .hap format spec

    TARGET refers to the ID of a variant or haplotype. LD is computed pair-wise between
    TARGET and all of the other haplotypes in the .hap (or genotype) file

    If TARGET is a variant ID, the ID must appear in GENOTYPES. Otherwise, it must
    be present in the .hap file
    """
    import logging

    from .ld import calc_ld

    log = logging.getLogger("haptools ld")
    db_time = "|%(asctime)s" if verbosity == "DEBUG" else ""
    logging.basicConfig(
        format="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
        datefmt="%H:%M:%S",
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

    if ids_file:
        with ids_file as id_file:
            ids = tuple(id_file.read().splitlines())
    elif ids:
        ids = tuple(ids)
    else:
        ids = None

    calc_ld(
        target,
        genotypes,
        haplotypes,
        region,
        samples,
        ids,
        chunk_size,
        discard_missing,
        from_gts,
        output,
        log,
    )


@main.command(short_help="Sort and index .hap files")
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--sort/--no-sort",
    is_flag=True,
    default=True,
    help="Sorting of the file will not be performed",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    show_default="input file",
    help="A .hap file containing sorted and indexed haplotypes and variants",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default="only errors",
    help="The level of verbosity desired",
)
def index(
    haplotypes: Path,
    sort: bool = False,
    output: Path = None,
    verbosity: str = "CRITICAL",
):
    """
    Takes in an unsorted .hap file and outputs it as a .gz and a .tbi file
    """

    import logging
    from .index import index_haps

    log = logging.getLogger("haptools index")
    logging.basicConfig(
        format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
    )

    index_haps(haplotypes, sort, output, log)


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")
