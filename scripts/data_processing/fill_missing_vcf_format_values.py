#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Fill missing FORMAT values in a VCF on a sample-specific basis
"""


import argparse
import pandas as pd
import pysam
from sys import stdin, stdout


type_map = {'Integer' : int, 'Float' : float, 'Character' : str, 'String' : str}


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')
    parser.add_argument('lookup_tsv', help='two-column .tsv of sample, value ' +
                        'indicating values to be filled')
    parser.add_argument('-f', '--field', default='GQ', help='Field to fill ' +
                        '[default: GQ]')
    parser.add_argument('vcf_out', help='output .vcf')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    # Determine type of variable for --field
    vmd = [v for k, v in invcf.header.formats.items() if k == args.field]
    if len(vmd) == 0:
        msg = 'Unable to find --field {} in input VCF header. Exiting.'
        exit(msg.format(args.field))
    else:
        vtype = vmd[0].type

    # Load lookup table as dict
    lt = pd.read_csv(args.lookup_tsv, sep='\t', header=None).set_index(0).to_dict()[1]

    # Open connection to output file
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=invcf.header)

    # Iterate over records in invcf, fill missing values, and write to vcf_out
    for record in invcf.fetch():
        for sid, val in lt.items():
            if record.samples[sid][args.field] is None:
                record.samples[sid][args.field] = type_map[vtype](val)
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

