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
from sys import stdin, stdout



def parse_vep_map(invcf):
    """
    Parse VEP field mappings to variable names
    """

    vep_text = invcf.header.info.get('CSQ').description
    vep_fields = vep_text.split('Format: ')[1].replace('\'', '').split('|')

    return {i : k for i, k in enumerate(vep_fields)}


def reformat_header(invcf):
    """
    Reformat header for output VCF
    """

    out_header = invcf.header

    return out_header


def cleanup(record, vep_map):
    """
    Simplify VEP entry for a single record
    """

    vep_vals = {}
    for i, vep_str in enumerate(record.info.get('CSQ')):
        vep_vals[i] = {k : v for k, v in zip(vep_map.values(), vep_str.split('|'))}

    

    import pdb; pdb.set_trace()



def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input VEP-annotated .vcf [default: stdin]',
                        default='stdin')
    parser.add_argument('vcf_out', help='output .vcf [default: stdout]', default='stdout')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Parse mapping of VEP fields
    vep_map = parse_vep_map(invcf)
    
    # Open connection to output file
    out_header = reformat_header(invcf)
    if args.vcf_out in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=out_header)

    # Iterate over records in invcf, clean up each, and write to outvcf
    for record in invcf.fetch():
        record = cleanup(record, vep_map)
        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

