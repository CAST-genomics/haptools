#!/usr/bin/env python

import sys
import pickle
import shutil
from pathlib import Path
from time import process_time

import click
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from haptools import logging
from haptools.data import GenotypesVCF, GenotypesPLINK


# COMMAND FOR GENERATING UKB PLOT:
# tests/bench_genotypes.py \
# --default-variants 18472 --default-samples 487409 --intervals-variants 1 80 4 \
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


def time_vcf(vcf, max_variants, chunk_size=500):
    GenotypesVCF(vcf).read(max_variants=max_variants)


def time_plink(pgen, max_variants, chunk_size=500):
    GenotypesPLINK(pgen).read(max_variants=max_variants)


def time_plink_chunk(pgen, max_variants, chunk_size=500):
    GenotypesPLINK(pgen, chunk_size=chunk_size).read(max_variants=max_variants)


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
    "-t",
    "--temp",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    default=Path("bench_temp_dir"),
    show_default=True,
    help="A temp directory into which to place generated temporary files",
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
    "-p",
    "--progress",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to display a progress bar. Useful for large benchmarks",
)
@click.option(
    "-s",
    "--silent",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to be entirely silent",
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
    temp,
    default_variants,
    default_samples,
    intervals_variants,
    intervals_samples,
    reps,
    output,
    progress=False,
    silent=False,
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
    LOG = logging.getLogger("run", 0 if silent else "DEBUG")
    gts = GenotypesPLINK(None, log=LOG)
    num_variants = max(DEFAULT_VARIANTS, INTERVALS_VARIANTS.stop)
    sample_size = max(DEFAULT_SAMPLES, INTERVALS_SAMPLES.stop)
    gts.samples = tuple(f"sample{i}" for i in range(sample_size))
    gts.variants = np.array(
        [(f"id{i}", "chr0", i, ("A", "T")) for i in range(1, num_variants + 1)],
        dtype=gts.variants.dtype,
    )
    np.random.seed(12345)
    if not temp.exists():
        LOG.info("Generating fake genotypes")
        gts.data = np.random.choice(
            [True, False], size=sample_size * num_variants * 2
        ).reshape((sample_size, num_variants, 2))
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
    FILE_TYPES = {
        "vcf": "VCF",
        "pgen": "PLINK2",
        "chunked200": "PLINK2 chunk_size: 200",
        "chunked400": "PLINK2 chunk_size: 400",
        "chunked600": "PLINK2 chunk_size: 600",
        "chunked800": "PLINK2 chunk_size: 800",
    }

    # create the files we will try to load if they haven't been created already
    if not gts.fname.exists():
        LOG.info("Creating VCF and PGEN files that we can load")
        variant_dir = create_variant_files(gts, INTERVALS_VARIANTS, DEFAULT_SAMPLES)
        sample_dir = create_sample_files(gts, INTERVALS_SAMPLES, DEFAULT_VARIANTS)
    else:
        variant_dir = gts.fname / "variant"
        sample_dir = gts.fname / "sample"
        LOG.info(
            "Temp directory already exists. Assuming files have already been created",
        )

    LOG.setLevel(level="ERROR")

    # run each test
    LOG.info("Benchmarking the loading of each file")
    results = {}
    len_variants = len(VARIANTS)
    for arg in ("samples", "variants"):
        genotype_dir = variant_dir
        if arg == "samples":
            genotype_dir = sample_dir
        results[arg] = {}
        intervals = INTERVALS_SAMPLES if arg == "samples" else INTERVALS_VARIANTS
        item_iter = progressbar(intervals, prefix=f"{arg}, {file_type}: ", out=sys.stderr) if progress else intervals
        for file_type in FILE_TYPES.keys():
            results[arg][file_type] = []
            for val in item_iter:
                chunk_size = 500
                if file_type == "vcf":
                    func = time_vcf
                    file = genotype_dir / f"{val}.vcf"
                elif file_type == "pgen" or file_type.startswith("chunked"):
                    func = time_plink if file_type == "pgen" else time_plink_chunk
                    file = genotype_dir / f"{val}.pgen"
                    if file_type.startswith("chunked"):
                        chunk_size = int(file_type[len("chunked") :])
                else:
                    continue
                times = np.empty(REPS, dtype=np.float64)
                for rep in range(REPS):
                    start = process_time()
                    func(file, max_variants=len_variants, chunk_size=chunk_size)
                    end = process_time()
                    times[rep] = end - start
                results[arg][file_type].append(times.mean())

    # plot the results
    LOG.info("Generating plot of results")
    fig, (ax_samples, ax_variants) = plt.subplots(1, 2, figsize=(10, 5))
    lab_font = {"fontsize": "xx-small"}
    for file_type in FILE_TYPES.keys():
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
    for file_type in FILE_TYPES.keys():
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
    fig.legend(loc="lower left", fontsize="xx-small")
    fig.tight_layout()
    fig.savefig(str(output), bbox_inches="tight", dpi=400)

    # saving results
    if archive is not None:
        with open(archive, "wb") as fh:
            pickle.dump(results, fh)


def test_bench_genotypes():
    tmp_dir = Path("bench_temp_dir")
    tmp_plot = Path("bench-plink.png")

    try:
        main(["-t", tmp_dir, "-o", tmp_plot, "--reps", 1, "-s"])
    except SystemExit as err:
        # re-raise unless main() finished without an error
        if err.code:
            raise

    shutil.rmtree(tmp_dir)
    tmp_plot.unlink()


if __name__ == "__main__":
    main()
