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

    prefixes = set(list(pop_map.values()) + list(cancer_map.values()))

    descriptions = {'AN' : 'Allele number',
                    'AC' : 'Allele count',
                    'AF' : 'Allele frequency',
                    'N_GT' : 'Number of informative genotypes',
                    'N_HOMREF' : 'Number of homozygous reference genotypes',
                    'N_HET' : 'Number of heterozygous/hemizygous genotypes',
                    'N_HOMALT' : 'Number of homozygous alternate genotypes',
                    'CARRIER_FREQ' : 'Fraction of individuals carrying non-reference genotypes'}
    float_keys = 'AF CARRIER_FREQ'.split()

    for prefix in [None] + list(prefixes):
        if prefix is None:
            prefix = ''
            descr_suffix = ''
        else:
            descr_suffix = ' in {} samples'.format(prefix)
            prefix = prefix + '_'
        for key, descr in descriptions.items():
            if key in float_keys:
                key_type = 'Float'
            else:
                key_type = 'Integer'
            if prefix + key in header.info.keys():
                header.info.remove_header(prefix + key)
            header.add_meta(key='INFO', 
                            items=[('ID', prefix + key), ('Number', 'A'), 
                                   ('Type', key_type),
                                   ('Description', descr + descr_suffix)])

    return header


def annotate_AF(record, pop_map, cancer_map):
    """
    Annotate AF by population and cancer type
    """

    # Prep counters
    counter_fields = 'AC AN N_GT N_HOMREF N_HET N_HOMALT'.split()
    counters = {f : 0 for f in counter_fields}
    for prefix in set(list(pop_map.values()) + list(cancer_map.values())):
        for f in counter_fields:
            counters['_'.join([prefix, f])] = 0

    for sid, sgt in record.samples.items():

        # Extract necessary sample info and skip if GT is missing
        GT = sgt['GT']
        pop = pop_map.get(sid)
        cancer = cancer_map.get(sid)
        if all([a is None for a in GT]):
            continue
        GT_called = tuple([a for a in GT if a is not None])
        AN = len(GT)
        AC = np.nansum(np.array([a for a in GT if a is not None]) > 0)

        # Update allele counters
        counters['AN'] += AN
        counters['AC'] += AC
        for prefix in [pop, cancer]:
            if prefix is not None:
                counters[prefix + '_AN'] += AN
                counters[prefix + '_AC'] += AC

        # Update GT counters
        is_ref = pd.Series([a == 0 for a in GT_called])
        if is_ref.all():
            gt_label = 'HOMREF'
        elif ~is_ref.all():
            gt_label = 'HOMALT'
        else:
            gt_label = 'HET'
        counters['N_GT'] += 1
        counters['N_' + gt_label] += 1
        for prefix in [pop, cancer]:
            counters[prefix + '_N_GT'] += 1
            counters[prefix + '_N_' + gt_label] += 1

    # Compute derived frequencies
    counters = {key : int(value) for key, value in counters.items()}
    if counters['AN'] > 0:
        counters['AF'] = counters['AC'] / counters['AN']
    if counters['N_GT'] > 0:
        n_nonref = counters['N_HET'] + counters['N_HOMALT']
        counters['CARRIER_FREQ'] = n_nonref / counters['N_GT']
    for prefix in set(list(pop_map.values()) + list(cancer_map.values())):
        if counters[prefix + '_AN'] > 0:
            counters[prefix + '_AF'] = counters[prefix + '_AC'] / counters[prefix + '_AN']
        if counters[prefix + '_N_GT'] > 0:
            n_nonref = counters[prefix + '_N_HET'] + counters[prefix + '_N_HOMALT']
            counters[prefix + '_CARRIER_FREQ'] = \
                n_nonref / counters[prefix + '_N_GT']

    # Write all values to record INFO
    for key, value in counters.items():
        record.info[key] = value

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

