#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Clean up verbose VEP output for RASMod VCFs
"""


import argparse
import pysam
from sys import stdout


def reformat_header(invcf):
    """
    Reformat header for output VCF
    """

    out_header = copy(invcf.header)

    return out_header


def cleanup(record):
    """
    Simplify VEP entry for a single record
    """

    import pdb; pdb.set_trace()



def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', help='input VEP-annotated .vcf', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='output .vcf ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Open connection to input vcf
    invcf = pysam.VariantFile(args.vcf)
    out_header = reformat_header(invcf)
    
    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.outfile, 'w', header=out_header)

    # Iterate over records in invcf, clean up each, and write to outvcf
    for record in invcf.fetch():
        record = cleanup(record)
        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

