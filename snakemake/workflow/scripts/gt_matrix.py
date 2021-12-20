#!/usr/bin/env python

import sys
import argparse
import pandas as pd


OTHER_COLS = ["ID", "MAF", "alleles"]


def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Use a genotype matrix to sum allele differences from the REF for STRs."
        )
    )
    parser.add_argument(
        "-o",
        "--out",
        default=sys.stdout,
        help="path to TSV file where variants are rows and samples are cols",
    )
    parser.add_argument(
        "-m",
        "--min-maf",
        type=restricted_float,
        default=0,
        help="discard all SNPs that have an MAF below this number",
    )
    parser.add_argument(
        "gt_matrix",
        nargs="?",
        default=sys.stdin,
        help=(
            "a tab-separated genotype matrix where variants are rows and samples are"
            " cols"
        ),
    )
    args = parser.parse_args()
    return args


def index_list_of_lists_by_list(list_of_lists, indices):
    for idx, items in zip(indices, list_of_lists):
        idx = idx[1]
        yield (items[idx[0]], items[idx[1]])


def create_snp_gt_matrix(snp_gt, min_maf):
    """convert the genotype entires into 0, 1, or 2 and threshold by MAF"""
    # which columns are the ones we need to change? all except the OTHER_COLS
    sample_cols = snp_gt.columns.difference(OTHER_COLS)
    # first, filter by MAF and replace the GTs with (0, 1, or 2)
    new_matrix = snp_gt[snp_gt["MAF"] > min_maf]
    new_matrix.loc[:, sample_cols] = new_matrix.loc[:, sample_cols].replace(
        ["0|0", "0|1", "1|0", "1|1"], [0, 1, 1, 2]
    )
    # now, we must flip entries where the alt allele is minor
    # we first retrieve the rows of the entries that need to be flipped by calculating
    # the freq of the alt allele and thresholding by 0.5, and then we just flip them
    alt_is_minor = (new_matrix[sample_cols].sum(axis=1) / (2 * len(sample_cols))) > 0.5
    new_matrix.loc[alt_is_minor, sample_cols] = new_matrix.loc[
        alt_is_minor, sample_cols
    ].replace({0: 2, 1: 1, 2: 0})
    return new_matrix


def create_str_gt_matrix(str_gt):
    """convert the genotype entries into genotype dosages"""
    str_gt.alleles = str_gt.alleles.str.split(",")
    # convert these to lengths
    alleles = [[len(al) for al in alls] for alls in str_gt.alleles.values.tolist()]
    # get the length of the ref allele, in each case
    ref_len = str_gt.alleles.str[0].str.len()
    for sample in str_gt.columns:
        if sample in OTHER_COLS:
            continue
        # note: this strategy is probably slow -- it creates a list of lists instead
        # of vectorizing
        samp_gt = str_gt[sample].str.split(r"\||\/", expand=True).astype("uint8")
        # convert to lengths
        samp_gt = pd.DataFrame(
            index_list_of_lists_by_list(alleles, samp_gt.iterrows()),
            index=ref_len.index,
        )
        # THIS IS WHERE THE MAGIC HAPPENS!
        # the GT is the sum of the differences between sample GT len and the REF len
        str_gt[sample] = samp_gt.sub(ref_len, axis=0).sum(axis=1)
    return str_gt


def main(args):
    gt = pd.read_csv(args.gt_matrix, sep="\t", index_col=0)

    # split the matrix into STRs and SNPs and drop SNPs that aren't bi-allelic
    gt["ID"] = gt["ID"].str.startswith("STR_")
    snp_gt = gt.loc[(~gt["ID"]) & (gt.alleles.str.len() > 2)]

    # convert to proper genotype values
    gt = pd.concat(
        [
            create_snp_gt_matrix(snp_gt, args.min_maf),
            create_str_gt_matrix(gt[gt["ID"]].copy()),
        ]
    ).sort_index()

    # convert ID column to snp_status and remove unnecessary alleles col
    gt.index = gt.index.astype(str) + ":" + gt["ID"].astype("uint8").astype(str)
    gt.drop(OTHER_COLS, axis=1, inplace=True)

    gt = gt.transpose()
    gt.index.rename("sample", inplace=True)
    gt.to_csv(args.out, sep="\t")


if __name__ == "__main__":
    main(parse_args())
