#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Annotate AFs by cancer type and population
"""


import argparse
import numpy as np
import pandas as pd
import pysam
from sys import stdin, stdout


def load_sample_metadata(meta_in, sample_col):
    """
    Load sample metadata necessary for AF annotation
    """

    df = pd.read_csv(meta_in, sep='\t')
    df.rename(columns={df.columns[0] : df.columns[0].replace('#', '')}, inplace=True)

    df.set_index(sample_col, drop=True, inplace=True)

    pop_map = df.POPULATION.to_dict()

    cancer_map = df.CANCER_TYPE.to_dict()

    return pop_map, cancer_map


def reformat_header(invcf, pop_map, cancer_map):
    """
    Add necessary columns to header
    """

    header = invcf.header

    return header


def annotate_AF(record, pop_map, cancer_map):
    """
    Annotate AF by population and cancer type
    """

    # Prep counters
    counter_fields = 'AC AN N_GT N_HOMREF N_HET N_HOMALT'.split()
    counters = {f : 0 for f in counter_fields}
    for pop in set(pop_map.values() + cancer_map.values()):
        for f in counter_fields:
            counters['_'.join([f, pop])] = 0

    for sid, sgt in record.samples.items():

        # Extract necessary sample info and skip if GT is missing
        GT = sgt['GT']
        pop = pop_map.get(sid)
        cancer = cancer_map.get(sid)
        if None in GT:
            continue
        import pdb; pdb.set_trace()
        AN = len(GT)
        AC = np.nansum(GT)

        # Update allele counters
        counters['AN'] += AN
        counters['AC'] += AC
        for suf in [pop, cancer]:
            if suf is not None:
                counters['AN_' + suf] += AN
                counters['AC' + suf] += AC

        # Update GT counters
        # TODO: implement this


    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf [default: stdin]', default='stdin')
    parser.add_argument('vcf_out', help='output .vcf [default: stdout]', default='stdout')
    parser.add_argument('-m', '--sample-metadata', required=True, help='.tsv ' + \
                        'of sample metadata. Must contain columns "POPULATION" ' + \
                        'and "CANCER_TYPE".')
    parser.add_argument('-s', '--sample-id-column', default='DONOR_ID',
                        help='Column to use for sample/participant ID [default: %default]')
    args = parser.parse_args()

    # Load sample population & cancer type info from --sample-metadata
    pop_map, cancer_map = \
        load_sample_metadata(args.sample_metadata, args.sample_id_column)

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Open connection to output file
    out_header = reformat_header(invcf, pop_map, cancer_map)
    if args.vcf_out in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=out_header)

    # Iterate over records in invcf, annotate AF, and write to outvcf
    for record in invcf.fetch():
        record = annotate_AF(record, pop_map, cancer_map)
        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

