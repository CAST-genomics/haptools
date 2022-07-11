#!/usr/bin/env python

import sys
import pickle
import logging
from pathlib import Path
from time import process_time

import click
import numpy as np
import matplotlib.pyplot as plt

from haptools.data import GenotypesRefAlt, GenotypesPLINK


# DEFAULT_SAMPLES = 500000
# # DEFAULT_VARIANTS = 20000
# # INTERVALS_VARIANTS = range(500, 20000, 500)
# DEFAULT_VARIANTS = 400
# INTERVALS_VARIANTS = range(10, 400, 10)
# INTERVALS_SAMPLES = range(12500, 500000, 12500)

DEFAULT_SAMPLES = 40
DEFAULT_VARIANTS = 40
INTERVALS_VARIANTS = range(5, 100, 5)
INTERVALS_SAMPLES = range(5, 100, 5)

# DEFAULT_SAMPLES = 5
# DEFAULT_VARIANTS = 4
# INTERVALS_VARIANTS = range(1, 4, 1)
# INTERVALS_SAMPLES = range(1, 5, 1)

REPS = 5
DATADIR = Path(__file__).parent.joinpath("data")


def create_variant_files(gts, intervals, num_samps):
    samples = gts.samples[:num_samps]
    variant_dir = gts.fname / "variant"
    variant_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        variants = tuple(gts.variants["id"][:val])
        sub = gts.subset(samples=samples, variants=variants)
        sub.fname = variant_dir / f"{val}.pgen"
        sub.write(clean_up=False)
    return variant_dir


def create_sample_files(gts, intervals, num_vars):
    variants = tuple(gts.variants["id"][:num_vars])
    sample_dir = gts.fname / "sample"
    sample_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        samples = gts.samples[:val]
        sub = gts.subset(samples=samples, variants=variants)
        sub.fname = sample_dir / f"{val}.pgen"
        sub.write(clean_up=False)
    return sample_dir


def time_vcf(vcf, max_variants):
    GenotypesRefAlt(vcf).read(max_variants=max_variants)


def time_plink(pgen, max_variants):
    GenotypesPLINK(pgen).read(max_variants=max_variants)


def time_plink_chunk(pgen, max_variants):
    GenotypesPLINK(pgen).read(max_variants=max_variants, chunk_size=50)


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
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
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
def main(pgen, temp, region, output, archive=None):
    """
    Benchmarks classes in the data.genotypes module

    GENOTYPES cab be formatted as a VCF or PLINK2 PGEN file

    \f
    Examples
    --------
    >>> haptools transform tests/data/example.vcf.gz tests/data/example.hap.gz > example_haps.vcf

    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF or PLINK2 PGEN format
    output : Path, optional
        The location to which to write output
    """
    global DEFAULT_VARIANTS, DEFAULT_SAMPLES, INTERVALS_VARIANTS, INTERVALS_SAMPLES
    print("Loading genotypes from PGEN file", file=sys.stderr)
    gts = GenotypesPLINK(pgen)
    gts.read(region=region, max_variants=DEFAULT_VARIANTS)
    gts.check_missing(discard_also=True)
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
    LOG = logging.getLogger("run")
    logging.basicConfig(
        format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
        level="ERROR",
    )

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
                    file = genotype_dir / f"{val}.vcf.gz"
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
    for file_type in ("vcf", "pgen", "chunked"):
        ax_samples.plot(
            INTERVALS_SAMPLES,
            results["samples"][file_type],
            marker="o",
            label=FILE_TYPES[file_type],
        )
    ax_samples.set_xlabel("Number of samples")
    ax_samples.set_ylim(ymin=0)
    for file_type in ("vcf", "pgen", "chunked"):
        ax_variants.plot(
            INTERVALS_VARIANTS,
            results["variants"][file_type],
            marker="o",
        )
    ax_variants.set_xlabel("Number of variants")
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
