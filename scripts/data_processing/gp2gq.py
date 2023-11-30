#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert genotype prior likelihoods (GP) to genotype qualities (GQ)
"""


import argparse
import pysam
from numpy import log10
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')
    parser.add_argument('vcf_out', help='output .vcf')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    # Open connection to output file
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=invcf.header)

    # Iterate over records in invcf
    for record in invcf.fetch():

        # Iterate over all samples per record
        for sid in samples:
            GT = record.samples[sid].get('GT', (None, None, ))
            AC = sum(GT)
            if not all([a is None for a in GT]):
                GP = record.samples[sid].get('GP', (None, None, None))[AC]
                if GP is None:
                    continue
                if GP == 1:
                    GQ = 99
                else:
                    GQ = int(round(-10 * log10(1 - GP), 0))
                record.samples[sid]['GQ'] = GQ

        # Write record to output VCF
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

