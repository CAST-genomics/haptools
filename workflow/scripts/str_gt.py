#!/usr/bin/env python

import sys
import argparse
import pandas as pd

from IPython import embed


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        "Use a genotype matrix to sum allele differences from the REF for STRs."
    )
    parser.add_argument(
        '-o', '--out', default=sys.stdout,
        help='path to TSV file where variants are rows and samples are cols'
    )
    parser.add_argument(
        'gt_matrix', nargs='?', default=sys.stdin,
        help='a tab-separated genotype matrix where variants are rows and samples are cols'
    )
    args = parser.parse_args()
    return args

def index_list_of_lists_by_list(list_of_lists, indices):
    for idx, items in zip(indices, list_of_lists):
        idx = idx[1]
        yield (items[idx[0]], items[idx[1]])


def main(args):
    gt = pd.read_csv(args.gt_matrix, sep="\t", index_col=0)
    # this is probably slow: it creates a list of lists instead of vectorizing
    gt.alleles = gt.alleles.str.split(',')
    # convert these to lengths
    alleles = [
        [len(al) for al in alls]
        for alls in gt.alleles.values.tolist()
    ]
    # get the length of the ref allele, in each case
    ref_len = gt.alleles.str[0].str.len()
    for sample in gt.columns:
        if sample == 'alleles':
            continue
        samp_gt = gt[sample].str.split(r'\||\/', expand=True).astype('uint8')
        # convert to lengths
        samp_gt = pd.DataFrame(index_list_of_lists_by_list(alleles, samp_gt.iterrows()), index=ref_len.index)
        # THIS IS WHERE THE MAGIC HAPPENS!
        # the GT is the sum of the differences between sample GT len and the REF len
        gt[sample] = samp_gt.sub(ref_len, axis=0).abs().sum(axis=1)
    gt.drop('alleles', axis=1, inplace=True)
    gt.to_csv(args.out, sep="\t")


if __name__ == '__main__':
    main(parse_args())
