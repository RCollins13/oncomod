#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Summarize number of somatic endpoints to test for each cancer type & gene
"""


import argparse
import os
import pandas as pd
import re
from sys import stdout, path
path.insert(0, os.path.join(path[0], '..', 'data_processing'))
from add_variant_set_members import load_members, get_members
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from general_utils import load_tx_map



# Define various variables used throughout the below functions
cancers = ['LUAD']
egfr_genes = ['EGFR']
egfr_chroms = {'7' : 'EGFR'}
category_descriptions = \
    {'mutations'        : '1. Frequent EGFR alterations',
     'codons'           : '2. Recurrently mutated EGFR codons',
     'exons'            : '3. Recurrently mutated EGFR exons',
     'burden'           : '4. Collapsed mutation sets',
     'comutations'      : '5. Frequent co-mutations involving EGFR',
     'egfr_nonegfr_comut' : '6. Frequent EGFR + non-EGFR comutations'}
tissue_map = {'LUAD' : 'lung'}


def infer_gene(set_id, tx_map):
    """
    Use a variety of strategies to infer EGFR gene relationship from set ID
    """

    if set_id.startswith('COMUT|'):
        genes = [infer_gene(g, tx_map) for g in set_id.split('|')[1:]]
        return ','.join(set([g for g in genes if g is not None]))
    elif set_id.startswith('ENST'):
        return tx_map['symbol'].get(set_id.split('_')[0])
    elif set_id.startswith('1_') \
         or set_id.startswith('11_') \
         or set_id.startswith('12_'):
        return egfr_chroms.get(set_id.split('_')[0])
    elif '_' in set_id:
        return set_id.split('_')[0]
    else:
        # If this point has been reached, need a new strategy to infer genes
        print('Unsure how to parse set ID ' + set_id + '\n')
        import pdb; pdb.set_trace()


def update_res(subres, infile, tx_map, members, combos_seen, min_freq=0.01):
    """
    Update category ID sets based on data in infile
    """

    # Load data from infile
    df = pd.read_csv(infile, sep='\t')

    # Attempt to infer gene for each category
    df['gene'] = df.set_id.apply(lambda x: infer_gene(x, tx_map))

    # Check for tissue-specific annotations, and set frequencies for non-matched
    # cancers to zero to skip the counting step
    for cancer, tissue in tissue_map.items():
        tissue_rows = df.set_id.str.contains(tissue)
        if tissue_rows.any():
            other_freq_cols = [c for c in df.columns if c.endswith('_AF') and cancer not in c]
            df.loc[tissue_rows, other_freq_cols] = 0

    # Deduplicate sets based on identical member variant IDs and frequencies
    df['vids'] = df.set_id.apply(get_members, members=members)
    df = df[~df.iloc[:, 1:].duplicated()]
    df.reset_index(inplace=True, drop=True)

    # Map categories onto cancer types & genes in subres
    for cancer in cancers:
        # Define a "hit" category per cancer type as one that exceeds min_freq
        # and hasn't already been seen in that same cancer type
        hits = (df[cancer + '_AF'] >= min_freq) & \
               ~(df.vids.isin(combos_seen[cancer]))
        combos_seen[cancer].update(set(df.vids[hits].tolist()))
        for idx, vals in df.loc[hits, 'set_id gene'.split()].iterrows():
            set_id, gstr = vals.values
            for gene in gstr.split(','):
                if gene not in subres[cancer].keys():
                    continue
                # Convert DEL/AMP to simple set ID to avoid slightly different variant IDs
                clean_parts = []
                for sub_id in set_id.split('|'):
                    if any([sub_id.endswith(x) for x in '_DEL _AMP _DUP'.split()]):
                        sub_id = '_'.join([gene, set_id.split('_')[-1]])
                    clean_parts.append(sub_id)
                if len(clean_parts) != len(set(clean_parts)):
                    continue
                clean_set_id = '|'.join(clean_parts)
                subres[cancer][gene].add(clean_set_id)

    # Return updated subres and updated list of member combos seen
    return subres, combos_seen


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutations', action='append', default=[],
                        help='frequencies corresponding to individual mutations')
    parser.add_argument('--codons', action='append', default=[],
                        help='frequencies corresponding to recurrent codon mutations')
    parser.add_argument('--exons', action='append', default=[],
                        help='frequencies corresponding to recurrent exon mutations')
    parser.add_argument('--burden-sets', action='append', default=[],
                        help='frequencies corresponding to burden sets')
    parser.add_argument('--comutations', action='append', default=[],
                        help='frequencies corresponding to comutation pairs')
    parser.add_argument('--egfr-nonegfr-comut', action='append', default=[],
                        help='frequencies corresponding to EGFR + non-EGFR ' + 
                        'comutation pairs')
    parser.add_argument('-t', '--transcript-info', required=True, help='.tsv ' + 
                        'mapping ENST:ENSG:symbol:length')
    parser.add_argument('--memberships', action='append', help='.tsv mapping ' +
                        'set IDs (first column) to constituent variant IDs ' +
                        '(final column). Can be specified multiple times.')
    parser.add_argument('-m', '--min-freq', default=0.01, type=float, help='Minimum ' + 
                        'frequency for a category to be retained per cancer type.')
    parser.add_argument('--min-egfr-nonegfr-comut', default=0.05, type=float, 
                        help='Minimum frequency for a EGFR-nonEGFR comutation pair ' +
                        'to be retained per cancer type.')
    parser.add_argument('-o', '--outfile', help='output .tsv of summary table ' + 
                        '[default: stdout]', default='stdout')
    parser.add_argument('-P', '--out-prefix', help='prefix for output lists of ' +
                        'somatic endpoints to test [default: %default]',
                        type=str, default='./')
    args = parser.parse_args()

    # Build dicts for collecting results
    res = {cat : {cncr : {gene : set() for gene in egfr_genes} for cncr in cancers} \
           for cat in category_descriptions.keys()}
    member_combos_seen = {cncr : set() for cncr in cancers}

    # Load mapping of set ID to constitutent variant IDs
    members = load_members(args.memberships)

    # Load transcript map
    tx_map = load_tx_map(args.transcript_info)
    
    # Load individual mutations
    for infile in args.mutations:
        res['mutations'], member_combos_seen = \
            update_res(res['mutations'], infile, tx_map, members, 
                       member_combos_seen, args.min_freq)

    # Load recurrently mutated codons
    for infile in args.codons:
        res['codons'], member_combos_seen = \
            update_res(res['codons'], infile, tx_map, members, 
                       member_combos_seen, args.min_freq)

    # Load recurrently mutated exons
    for infile in args.exons:
        res['exons'], member_combos_seen = \
            update_res(res['exons'], infile, tx_map, members, 
                       member_combos_seen, args.min_freq)

    # Load burden sets
    for infile in args.burden_sets:
        res['burden'], member_combos_seen = \
            update_res(res['burden'], infile, tx_map, members, 
                       member_combos_seen, args.min_freq)

    # Load comutation pairs
    for infile in args.comutations:
        res['comutations'], member_combos_seen = \
            update_res(res['comutations'], infile, tx_map, members, 
                       member_combos_seen, args.min_freq)

    # Load EGFR + non-EGFR comutation pairs
    for infile in args.egfr_nonegfr_comut:
        res['egfr_nonegfr_comut'], member_combos_seen = \
            update_res(res['egfr_nonegfr_comut'], infile, tx_map, members, 
                       member_combos_seen, args.min_egfr_nonegfr_comut)

    # Output lists of somatic endpoints per gene & cancer type
    for cancer in cancers:
        for gene in egfr_genes:
            fout = open('{}{}.{}.somatic_endpoints.tsv'.format(args.out_prefix, cancer, gene), 'w')
            for cat in category_descriptions.keys():
                for val in res[cat][cancer][gene]:
                    mems = get_members(val, members)
                    fout.write('{}\t{}\n'.format(val, mems))
            fout.close()

    # Open connection to --outfile
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    header_vals = ['Somatic Criteria']
    for cancer in cancers:
        for gene in egfr_genes:
            header_vals.append('_'.join([cancer, gene]))
        header_vals.append(cancer + '_Union')
    outfile.write('\t'.join(header_vals + ['Total']) + '\n')

    # Summarize as table
    for category, vals in res.items():
        outvals = [category_descriptions[category]]
        total = 0
        for cancer in cancers:
            sub_union = set()
            for gene in egfr_genes:
                outvals.append(len(vals[cancer][gene]))
                sub_union.update(vals[cancer][gene])
            outvals.append(len(sub_union))
            total += len(sub_union)
        outfile.write('\t'.join([str(x) for x in outvals + [total]]) + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

