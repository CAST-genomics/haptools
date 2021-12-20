#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot phenotypes against genotypes in a phens.tsv.gz file."
    )
    parser.add_argument("-o", "--out", help="path to plot PNG file")
    parser.add_argument(
        "-i",
        "--index",
        default=0,
        type=int,
        help="the index of the column in the matrix that contains the genotypes",
    )
    parser.add_argument(
        "phens",
        nargs="?",
        default=sys.stdin,
        help=(
            "a tab-separated table composed of genotypes and a phenotype (the last"
            " column)"
        ),
    )
    args = parser.parse_args()
    return args


def main(args):
    phens = pd.read_csv(args.phens, sep="\t", index_col=0)
    phen_plot = sns.regplot(x=phens.iloc[:, args.index], y=phens.iloc[:, -1])
    phen_plot.get_figure().savefig(args.out)


if __name__ == "__main__":
    main(parse_args())
