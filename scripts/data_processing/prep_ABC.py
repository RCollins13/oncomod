#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Curate ABC promoters for VEP annotation
"""


import argparse
import pandas as pd
from sys import stdout


def get_abc_iterables(tsv_in):
    """
    Load & filter ABC enhancers
    Returns a generator of pd.Series, one per enhancer-gene pair
    """

    # Load & clean data
    keep_cols = 'chr start end TargetGene ABC.Score CellType'
    abc = pd.read_csv(tsv_in, sep='\t', usecols=keep_cols.split())
    abc.chr = abc.chr.str.replace('chr', '')

    # Subset to cell types of interest
    tissue_map = {'pancreas-Roadmap' : 'pancreas',
                  'Panc1-ENCODE' : 'pancreas',
                  'body_of_pancreas-ENCODE' : 'pancreas',
                  'sigmoid_colon-ENCODE' : 'colon',
                  'transverse_colon-ENCODE' : 'colon',
                  'fibroblast_of_lung-Roadmap' : 'lung'}
    abc = abc[abc.CellType.isin(tissue_map.keys())]
    abc['Tissue'] = abc.CellType.map(tissue_map)

    # Return as generator of row Series
    return abc.iterrows()


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--abc', help='Input ABC .tsv', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='output .bed ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Iterate over features in input file
    ofmt = '{}\t{}\t{}\t{}:{}:{}:{}\n'
    for idx, edata in get_abc_iterables(args.abc):
        outfile.write(ofmt.format(*edata.values))

    # Close connection to outfile to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

