#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert a VCF into a simple allele dosage matrix
"""


import argparse
import numpy as np
import pandas as pd
import pysam
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')
    parser.add_argument('outfile', help='output .tsv')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('\t'.join(['#VID'] + samples) + '\n')

    # Iterate over records in invcf, convert to dosage vectors, and write to outfile
    for record in invcf.fetch():
        outvals = [record.id]
        for sid, sgt in record.samples.items():
            GT = [a for a in sgt['GT'] if a is not None]
            if len(GT) == 0:
                outvals.append('NA')
            else:
                dos = int(np.sum([a for a in GT if a > 0]))
                outvals.append(str(dos))
        outfile.write('\t'.join(outvals) + '\n')

    # Close connection to output file to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

