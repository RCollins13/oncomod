#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Filter pre-computed variant/set frequency tables
"""


import argparse
import pandas as pd
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--freq-tsv', default='stdin', help='variant frequency ' + \
                        '.tsv [default: stdin]')
    parser.add_argument('-m', '--min-freq', default=0, type=float, help='Minimum ' + \
                        'frequency in any cancer type to be retained.')
    parser.add_argument('--min-ac', default=0, type=int, help='Minimum allele ' + \
                        'count in any cancer type to be retained.')
    parser.add_argument('--report-ac', default=False, action='store_true',
                        help='report AC values per row in --outfile [default: report AF]')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Load frequency table
    if args.freq_tsv in '- stdin /dev/stdin'.split():
        fin = stdin
    else:
        fin = args.freq_tsv
    df = pd.read_csv(fin, sep='\t')

    # Filter on any frequency >= --min-freq and AC >= --min-ac
    af_cols = [c for c in df.columns if c.endswith('_AF')]
    df['max_freq'] = df.loc[:, af_cols].max(axis=1)
    ac_cols = [c for c in df.columns if c.endswith('_AC')]
    df['max_AC'] = df.loc[:, ac_cols].max(axis=1)
    df.sort_values('max_freq', inplace=True, ascending=False)
    keep = (df.max_freq >= args.min_freq) & (df.max_AC >= args.min_ac)
    if args.report_ac:
        out_df = df.loc[keep, ['set_id'] + ac_cols]
    else:
        out_df = df.loc[keep, ['set_id'] + af_cols]

    # Write to --outfile
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = args.outfile
    out_df.to_csv(outfile, header=True, index=False, sep='\t')


if __name__ == '__main__':
    main()

