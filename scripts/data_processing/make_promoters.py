#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Create a BED file of gene promoters from an input GTF
"""


import argparse
import pybedtools as pbt
from sys import stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', help='Input GTF', required=True)
    parser.add_argument('--distance', default=2000, type=int, help='Distance ' +
                        'upstream of TSS to define as a promoter')
    parser.add_argument('--coding-only', default=False, action='store_true',
                        help='Restrict to protein-coding genes')
    parser.add_argument('-o', '--outfile', default='stdout', help='output .bed ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Iterate over gene features in gtf
    for feature in pbt.BedTool(args.gtf).filter(lambda x: x[2] == 'gene'):

        if args.coding_only and feature.attrs.get('gene_type', 'unk') != 'protein_coding':
            continue
        
        chrom = str(feature.chrom)
        gene = feature.attrs.get('gene_name')
        tx = feature.attrs.get('transcript_id')
        
        if feature.strand == '+':
            end = int(feature.start)
            start = max([0, end - args.distance])
        elif feature.strand == '-':
            start = feature.end
            end = start + args.distance
        
        if not any(x is None for x in [gene, tx]):
            outfile.write('{}\t{}\t{}\t{}|{}\n'.format(chrom, start, end, gene, tx))

    # Close connection to outfile to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

