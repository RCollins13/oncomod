#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Unify TCGA sample/donors lists between exomes, genotyped arrays, and imputed arrays
"""


import argparse
import pandas as pd
from sys import stdout


def _parse_donor(id_str, skip=[]):
    """
    Extract TCGA donor ID from various TCGA-styled id_str inputs
    """

    return '-'.join([x for i, x in enumerate(id_str.split('-')) if i not in skip][:3])


def load_exome_samples(exome_ids_in, exome_id_map_in):
    """
    Load list of exome samples and split to donor-indexed dataframe
    """

    # Read list of exome IDs
    ex_df = pd.read_csv(exome_ids_in, sep='\t', header=None, names=['WES_BAM_ID'])

    # Read ID map
    id_map = pd.read_csv(exome_id_map_in, sep='\t', usecols='#WES_ID TCGA_sample_ID'.split())
    id_map.rename(columns={'#WES_ID' : 'WES_BAM_ID',
                          'TCGA_sample_ID' : 'WES_SAMPLE_ID'},
                  inplace=True)
    id_map['DONOR_ID'] = id_map.WES_SAMPLE_ID.apply(_parse_donor)
    id_map.drop_duplicates('DONOR_ID', keep='first', inplace=True)

    # Merge exome IDs and ID map
    return ex_df.merge(id_map, on='WES_BAM_ID', how='inner')


def load_array_samples(typed_in, imputed_in):
    """
    Load list of array samples and merge typed & imputed lists
    """

    # Load & donor-deduplicate typed input
    typed_df = pd.read_csv(typed_in, sep='\t', header=None, names=['ARRAY_TYPED_ID'])
    typed_df['DONOR_ID'] = typed_df.ARRAY_TYPED_ID.apply(_parse_donor)
    typed_df.drop_duplicates('DONOR_ID', keep='first', inplace=True)

    # Load & donor-deduplicated imputed input
    imputed_df = pd.read_csv(imputed_in, sep='\t', header=None, names=['ARRAY_IMPUTED_ID'])
    imputed_df['DONOR_ID'] = imputed_df.ARRAY_IMPUTED_ID.apply(_parse_donor, skip=[1])
    imputed_df.drop_duplicates('DONOR_ID', keep='first', inplace=True)

    # Merge typed & imputed inputs
    return typed_df.merge(imputed_df, on='DONOR_ID', how='left')


def format_outputs(ex_df, ar_df, tss_table_in, out_prefix):
    """
    Merges exome & array sample tables
    Infers cancer type for each donor
    Writes various output files
    """

    # Merge and reorder exome & array data
    col_order = 'DONOR_ID WES_SAMPLE_ID WES_BAM_ID ARRAY_TYPED_ID ARRAY_IMPUTED_ID'
    out_df = ex_df.merge(ar_df, on='DONOR_ID', how='outer')[col_order.split()]

    import pdb; pdb.set_trace()


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exome-ids', help='List of exome sample', required=True)
    parser.add_argument('--exome-id-map', help='.tsv mapping exome CRAM/BAM IDs ' +
                        'to their TCGA sample IDs', required=True)
    parser.add_argument('--array-typed-ids', help='List of array-genotyped ' + 
                        'samples', required=True)
    parser.add_argument('--array-imputed-ids', help='List of array-imputed ' + 
                        'samples', required=True)
    parser.add_argument('--tcga-tss-table', help='TCGA TSS code table', required=True)
    parser.add_argument('--out-prefix', help='Prefix for all output files', 
                        required=True)
    args = parser.parse_args()

    # Load exome samples
    ex_df = load_exome_samples(args.exome_ids, args.exome_id_map)

    # Load array samples
    ar_df = load_array_samples(args.array_typed_ids, args.array_imputed_ids)

    # Merge exome & array data and write to output files
    format_outputs(ex_df, ar_df, args.tcga_tss_table, args.out_prefix)



if __name__ == '__main__':
    main()

