#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Identify individual exons recurrently mutated at two or more distinct codons
"""


import argparse
import pandas as pd
from sys import stdin, stdout


set_id_fmt = '{}_exon_{}'


def _process_groupby(gbo):
    """
    Count number of unique codons in a pd.DataFrameGroupBy and concatenate their VIDs
    """

    vids = gbo.vids.str.split(',').values.tolist()
    vids = ','.join(set([vid for sub in vids for vid in sub]))
    codons = set(gbo.codon)

    return len(codons), vids


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('coding_tsv', help='collapsed coding consequences .tsv')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Load coding consequences tsv
    cdf = pd.read_csv(args.coding_tsv, sep='\t')
    cdf.exon = cdf.exon.fillna(0).astype(int)

    # Open connection to outfile
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('set_id\ttranscript\texon\tvids\n')
    # Group entries by transcript and exon
    group_res = cdf.groupby(by='transcript exon'.split()).apply(_process_groupby)
    for info, res in group_res.iteritems():
        if res[0] > 1:
            outvals = [set_id_fmt.format(*info), *info, res[1]]
            outfile.write('\t'.join([str(x) for x in outvals]) + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

