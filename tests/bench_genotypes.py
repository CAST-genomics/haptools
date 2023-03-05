#!/usr/bin/env python

import sys
import pickle
import logging
from pathlib import Path
from time import process_time

import click
import numpy as np
import matplotlib.pyplot as plt

from haptools.data import GenotypesVCF, GenotypesPLINK


# COMMAND FOR GENERATING UKB PLOT:
# tests/bench_genotypes.py -p /projects/ps-gymreklab/CAST/amassara/plink.pgen \
# --default-variants 20000 --default-samples 500000 --intervals-variants 1 80 4 \
# --intervals-samples 1 80 4 -o plot.png -a results.pickle

DATADIR = Path(__file__).parent.joinpath("data")


def create_variant_files(gts, intervals, num_samps):
    samples = gts.samples[:num_samps]
    variant_dir = gts.fname / "variant"
    variant_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        variants = tuple(gts.variants["id"][:val])
        sub = gts.subset(samples=samples, variants=variants)
        # write PLINK2 files
        sub.fname = variant_dir / f"{val}.pgen"
        sub.write()
        # write VCF files
        vcf = GenotypesVCF(sub.fname.with_suffix(".vcf"))
        vcf.data = sub.data
        vcf.samples = sub.samples
        vcf.variants = sub.variants
        vcf.write()
    return variant_dir


def create_sample_files(gts, intervals, num_vars):
    variants = tuple(gts.variants["id"][:num_vars])
    sample_dir = gts.fname / "sample"
    sample_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        samples = gts.samples[:val]
        sub = gts.subset(samples=samples, variants=variants)
        # write PLINK2 files
        sub.fname = sample_dir / f"{val}.pgen"
        sub.write()
        # write VCF files
        vcf = GenotypesVCF(sub.fname.with_suffix(".vcf"))
        vcf.data = sub.data
        vcf.samples = sub.samples
        vcf.variants = sub.variants
        vcf.write()
    return sample_dir


def time_vcf(vcf, max_variants):
    GenotypesVCF(vcf).read(max_variants=max_variants)


def time_plink(pgen, max_variants):
    GenotypesPLINK(pgen).read(max_variants=max_variants)


def time_plink_chunk(pgen, max_variants):
    GenotypesPLINK(pgen, chunk_size=50).read(max_variants=max_variants)


def progressbar(it, prefix="", size=60, out=sys.stdout):  # Python3.6+
    count = len(it)

    def show(j):
        x = int(size * j / count)
        print(
            f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}",
            end="\r",
            file=out,
            flush=True,
        )

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("\n", flush=True, file=out)


@click.command()
@click.option(
    "-p",
    "--pgen",
    type=click.Path(exists=True, path_type=Path),
    default=DATADIR.joinpath("simple.pgen"),
    show_default=True,
    help="A PGEN file containing the genotypes",
)
@click.option(
    "-t",
    "--temp",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    default=Path("bench_temp_dir"),
    show_default=True,
    help="A temp directory into which to place generated temporary files",
)
@click.option(
    "-r",
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help=(
        "The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'"
        "For this to work, the VCF must be indexed and the seqname provided must "
        "correspond with one in the files"
    ),
)
@click.option(
    "--default-variants",
    type=int,
    default=4,
    show_default=True,
    help="The number of variants to use when we vary the number of samples",
)
@click.option(
    "--default-samples",
    type=int,
    default=5,
    show_default=True,
    help="The number of samples to use when we vary the number of variants",
)
@click.option(
    "--intervals-variants",
    type=int,
    nargs=3,
    default=(1, 4, 1),
    show_default=True,
    help="The start, end, and step values for the x-axis of the variants plot",
)
@click.option(
    "--intervals-samples",
    type=int,
    nargs=3,
    default=(1, 5, 1),
    show_default=True,
    help="The start, end, and step values for the x-axis of the samples plot",
)
@click.option(
    "--reps",
    type=int,
    default=3,
    show_default=True,
    help="For each benchmark value, we take the mean of --reps X replicates",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("bench-plink.png"),
    show_default=True,
    help="A PNG file containing the results of the benchmark",
)
@click.option(
    "-a",
    "--archive",
    type=click.Path(path_type=Path),
    default=None,
    show_default="do not save generated results",
    help="A python pickle file into which to store results",
)
def main(
    pgen,
    temp,
    region,
    default_variants,
    default_samples,
    intervals_variants,
    intervals_samples,
    reps,
    output,
    archive=None,
):
    """
    Benchmarks classes in the data.genotypes module
    """
    DEFAULT_VARIANTS, DEFAULT_SAMPLES, INTERVALS_VARIANTS, INTERVALS_SAMPLES, REPS = (
        default_variants,
        default_samples,
        range(*intervals_variants),
        range(*intervals_samples),
        reps,
    )
    LOG = logging.getLogger("run")
    logging.basicConfig(
        format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
        level="DEBUG",
    )
    gts = GenotypesPLINK(pgen, chunk_size=500, log=LOG)
    if not temp.exists():
        print("Loading genotypes from PGEN file", file=sys.stderr)
        gts.read(
            region=region, max_variants=max(DEFAULT_VARIANTS, INTERVALS_VARIANTS.stop)
        )
        gts.check_missing(discard_also=True)
    else:
        gts.read_variants(max_variants=max(DEFAULT_VARIANTS, INTERVALS_VARIANTS.stop))
        gts.read_samples()
    gts.fname = temp
    # set initial variables
    SAMPLES = gts.samples
    VARIANTS = gts.variants["id"]
    DEFAULT_SAMPLES = min(DEFAULT_SAMPLES, len(SAMPLES))
    DEFAULT_VARIANTS = min(DEFAULT_VARIANTS, len(VARIANTS))
    if INTERVALS_VARIANTS.stop > len(VARIANTS):
        INTERVALS_VARIANTS = range(
            INTERVALS_VARIANTS.start, len(VARIANTS), INTERVALS_VARIANTS.step
        )
    INTERVALS_VARIANTS = list(INTERVALS_VARIANTS)
    if INTERVALS_SAMPLES.stop > len(SAMPLES):
        INTERVALS_SAMPLES = range(
            INTERVALS_SAMPLES.start, len(SAMPLES), INTERVALS_SAMPLES.step
        )
    INTERVALS_SAMPLES = list(INTERVALS_SAMPLES)
    FILE_TYPES = {"vcf": "VCF", "pgen": "PLINK2", "chunked": "PLINK2 chunked"}

    # create the files we will try to load if they haven't been created already
    if not gts.fname.exists():
        print("Creating VCF and PGEN files that we can load", file=sys.stderr)
        variant_dir = create_variant_files(gts, INTERVALS_VARIANTS, DEFAULT_SAMPLES)
        sample_dir = create_sample_files(gts, INTERVALS_SAMPLES, DEFAULT_VARIANTS)
    else:
        variant_dir = gts.fname / "variant"
        sample_dir = gts.fname / "sample"
        print(
            "Temp directory already exists. Assuming files have already been created",
            file=sys.stderr,
        )

    logging.getLogger().setLevel(level="ERROR")

    # run each test
    print("Benchmarking the loading of each file", file=sys.stderr)
    results = {}
    for arg in ("samples", "variants"):
        genotype_dir = variant_dir
        if arg == "samples":
            genotype_dir = sample_dir
        results[arg] = {}
        intervals = INTERVALS_SAMPLES if arg == "samples" else INTERVALS_VARIANTS
        for file_type in ("vcf", "pgen", "chunked"):
            results[arg][file_type] = []
            for val in progressbar(
                intervals, prefix=f"{arg}, {file_type}: ", out=sys.stderr
            ):
                if file_type == "vcf":
                    func = time_vcf
                    file = genotype_dir / f"{val}.vcf"
                elif file_type == "pgen" or file_type == "chunked":
                    func = time_plink if file_type == "pgen" else time_plink_chunk
                    file = genotype_dir / f"{val}.pgen"
                else:
                    continue
                times = np.empty(REPS, dtype=np.float64)
                for rep in range(REPS):
                    start = process_time()
                    func(file, max_variants=len(VARIANTS))
                    end = process_time()
                    times[rep] = end - start
                results[arg][file_type].append(times.mean())

    # plot the results
    print("Generating plot of results", file=sys.stderr)
    fig, (ax_samples, ax_variants) = plt.subplots(1, 2, figsize=(10, 5))
    lab_font = {"fontsize": "xx-small"}
    for file_type in ("vcf", "pgen", "chunked"):
        x_vals = INTERVALS_SAMPLES
        y_vals = results["samples"][file_type]
        # fit a line to each so that we can report the slope
        slope = np.polyfit(x_vals, y_vals, 1)[0]
        ax_samples.plot(
            x_vals,
            y_vals,
            marker="o",
            label=FILE_TYPES[file_type],
        )
        ax_samples.text(
            x_vals[-1],
            y_vals[-1] + (y_vals[-1] / 16),
            f"m = {slope:.3E}",
            fontdict=lab_font,
        )
    ax_samples.set_xlabel(f"Number of samples\nnum_variants = {DEFAULT_VARIANTS}")
    ax_samples.set_ylim(ymin=0)
    for file_type in ("vcf", "pgen", "chunked"):
        x_vals = INTERVALS_VARIANTS
        y_vals = results["variants"][file_type]
        # fit a line to each so that we can report the slope
        slope = np.polyfit(x_vals, y_vals, 1)[0]
        ax_variants.plot(
            x_vals,
            y_vals,
            marker="o",
        )
        ax_variants.text(
            x_vals[-1],
            y_vals[-1] + (y_vals[-1] / 16),
            f"m = {slope:.3E}",
            fontdict=lab_font,
        )
    ax_variants.set_xlabel(f"Number of variants\nnum_samples = {DEFAULT_SAMPLES}")
    ax_variants.set_ylim(ymin=0)
    fig.supylabel("CPU Time (s)")
    fig.legend(loc="lower left", fontsize="x-small")
    fig.tight_layout()
    fig.savefig(str(output), bbox_inches="tight", dpi=400)

    # saving results
    if archive is not None:
        with open(archive, "wb") as fh:
            pickle.dump(results, fh)


if __name__ == "__main__":
    main()
