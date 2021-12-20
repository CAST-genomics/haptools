#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import transforms
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Use a genotype matrix to create an LD heatmap."
    )
    parser.add_argument(
        "-o",
        "--out",
        default=sys.stdout,
        help="path to a PNG file to which to write the LD heatmap",
    )
    parser.add_argument(
        "gt_matrix",
        nargs="?",
        default=sys.stdin,
        help="a tab-separated GT matrix where variants are cols and samples are rows",
    )
    args = parser.parse_args()
    return args


def main(args):
    # import the data
    gt = pd.read_csv(args.gt_matrix, sep="\t", index_col=0)
    # calculate the correlation (r^2) between all columns (ie SNPs)
    corr = gt.corr().pow(2)
    # generate a mask for the lower diagonal (to create a triangle)
    mask = np.triu(np.ones_like(corr, dtype=bool))
    # create a 45 deg transformation
    # first of all, the base transformation of the data points is needed
    # base = plt.gca().transData
    # rot = transforms.Affine2D().rotate_deg(45)
    # create the heatmap and write it to the output file
    sns.heatmap(
        corr,
        mask=mask,
        xticklabels=False,
        yticklabels=False,
        cmap="YlGnBu",
        # transform=(rot+base)
    )
    plt.savefig(args.out)


if __name__ == "__main__":
    main(parse_args())
