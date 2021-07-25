#!/usr/bin/env python

import sys
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        "Use a genotype matrix to create phenotypes from a variant or variants."
    )
    parser.add_argument(
        '-o', '--out', default=sys.stdout,
        help='path to TSV file where variants are cols and samples are rows'
    )
    parser.add_argument(
        "--beta0", type=float, default=0,
        help="beta values in order from smallest to largest (ie beta0, beta1, beta2)"
    )
    parser.add_argument(
        "--beta1", nargs='+', type=float, default=[],
        help="beta values in order from smallest to largest (ie beta0, beta1, beta2)"
    )
    parser.add_argument(
        '--index', nargs='*',
        help='the index of the variant(s) to use when doing the simulations'
    )
    parser.add_argument(
        '--max-vars', type=int, help='the max number of variants to consider'
    )
    parser.add_argument(
        'gt_matrix', nargs='?', default=sys.stdin,
        help='a tab-separated genotype matrix where variants are cols and samples are rows'
    )
    args = parser.parse_args()
    return args


def main(args):
    if args.index:
        gt = pd.read_csv(args.gt_matrix, sep="\t", usecols=args.index, index_col=0)
    else:
        gt = pd.read_csv(args.gt_matrix, sep="\t")
        if args.max_vars:
            gt = gt.sample(args.max_vars)
    phen = args.betas[0]
    for beta_idx, beta in enumerate(args.betas[1:]):
        phen = ((beta**beta_idx) * gt.iloc[:,0]) + phen
    gt['phen'].to_csv(args.out, header=False, sep="\t")


if __name__ == '__main__':
    main(parse_args())
