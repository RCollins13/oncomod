#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Define well-captured exome intervals for a subset of samples of interest
"""


import argparse
import pandas as pd
import re
from sys import stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--coverage-matrix', help='Coverage matrix .tsv', required=True)
    parser.add_argument('--samples-list', help='List of samples to consider. ' +
                        'If not provided, will use all samples.')
    parser.add_argument('--min-frac-samples', help='Minimum fraction of samples ' +
                        'that must be covered at --min-frac-target. [default: %default]',
                        default=0.9, type=float)
    parser.add_argument('--min-frac-target', help='Minimum fraction of target ' +
                        'that must be covered in at least --min-frac-samples. ' +
                        '[default: %default]', default=0.9, type=float)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to output ' +
                        'BED file of well-captured intervals [default: stdout]')
    args = parser.parse_args()

    # Read list of samples of interest, if optioned
    if args.samples_list is not None:
        with open(args.samples_list) as fin:
            samples = set([s.rstrip() for s in fin.readlines()])
        usecols = list(samples)
    else:
        usecols = None

    # Load coverage matrix with intervals as index
    df = pd.read_csv(args.coverage_matrix, sep='\t', usecols=usecols,
                     index_col='Target', dtype=float)

    # Compute fraction of samples covered at --min-frac-target
    def _calc_frac(values, thresh=0.9):
        return (values >= thresh).sum() / len(values)
    frac_samples = df.apply(_calc_frac, thresh=args.min_frac_target, axis=1)

    # Identify targets covered in --min-frac-samples
    keep_targets = (frac_samples >= args.min_frac_samples)

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')

    # Split eligible targets and write to BED output file
    def _format_bed(istr):
        return '\t'.join(re.split(':|-', istr)) + '\n'
    for interval in df.index[keep_targets].apply(_format_bed):
        fout.write(interval)
    fout.close()


if __name__ == '__main__':
    main()

