#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pandas as pd
from math import sqrt
from numpy.random import normal


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        "Use a genotype matrix to create phenotypes from a variant or variants."
    )
    parser.add_argument(
        '-o', '--out', default=sys.stdout,
        help='path to TSV file where variants are cols, last column is phenotypes, and samples are rows'
    )
    mut_ex1 = parser.add_mutually_exclusive_group()
    mut_ex1.add_argument(
        "--beta-snp", type=float, default=0,
        help="beta value for a randomly chosen SNP"
    )
    mut_ex1.add_argument(
        "--beta-str", type=float, default=0,
        help="beta value for a randomly chosen STR"
    )
    mut_ex2 = parser.add_mutually_exclusive_group()
    mut_ex2.add_argument(
        '--str-loc', nargs='+',
        help='the start POS of the STR(s) to use for simulating the phenotypes'
    )
    mut_ex2.add_argument(
        '--snp-loc', nargs='+',
        help='the start POS of the SNP(s) to use for simulating the phenotypes'
    )
    mut_ex2.add_argument(
        '--max-vars', type=int, default=1, help='the max number of random variants to consider'
    )
    parser.add_argument(
        'gt_matrix', nargs='?', default=sys.stdin,
        help='a tab-separated GT matrix where variants are cols and samples are rows'
    )
    args = parser.parse_args()
    return args


def main(args):
    np.random.seed(40)
    if args.str_loc or args.snp_loc:
        variants = ['sample']
        if args.str_loc:
            variants.extend([
                idx+":1" for idx in args.str_loc
            ])
        else:
            variants.extend([
                idx+":0" for idx in args.snp_loc
            ])
        gt = pd.read_csv(
            args.gt_matrix, sep="\t", index_col=0, usecols=variants
        )
    elif args.max_vars:
        gt = pd.read_csv(
            args.gt_matrix, sep="\t", index_col=0
        )
        gt = gt.sample(args.max_vars, axis=1)
    else:
        # this shouldn't happen!
        pass

    gt['phen'] = 0
    for col in gt.columns:
        if col == 'phen':
            continue
        # z-normalize the column so it has stdev 1
        gt[col] = (gt[col] - gt[col].mean())/gt[col].std(ddof=0)
        # use the STR beta if the col name has a '1' at the end of it
        beta_val = args.beta_str if int(col[-1]) else args.beta_snp
        gt['phen'] = gt[col]*beta_val + gt['phen']
    # add some noise! sample randomly from a gaussian distribution
    gt['phen'] = normal(scale=sqrt(1-(beta_val**2)), size=gt[col].shape)

    gt.to_csv(args.out, header=True, sep="\t")


if __name__ == '__main__':
    main(parse_args())
