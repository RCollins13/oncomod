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
from sys import stdout, path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from general_utils import load_tx_map


# Define various variables used throughout the below functions
cancers = 'PDAC CRAD SKCM LUAD'.split()
ras_genes = 'NRAS HRAS KRAS'.split()
ras_chroms = {'1' : 'NRAS', '11' : 'HRAS', '12' : 'KRAS'}
category_descriptions = \
    {'mutations'        : '1. Frequent RAS alterations',
     'codons'           : '2. Recurrently mutated RAS codons',
     'burden'           : '3. Collapsed mutation sets',
     'comutations'      : '4. Frequent co-mutations involving RAS',
     'ras_nonras_comut' : '5. Frequent RAS + non-RAS comutations'}
tissue_map = {'PDAC' : 'pancreas',
              'CRAD' : 'colon',
              'SKCM' : 'skin',
              'LUAD' : 'lung'}


def infer_gene(set_id, tx_map):
    """
    Use a variety of strategies to infer RAS gene relationship from set ID
    """

    if set_id.startswith('COMUT|'):
        genes = [infer_gene(g, tx_map) for g in set_id.split('|')[1:]]
        return ','.join(set([g for g in genes if g is not None]))
    elif set_id.startswith('ENST'):
        return tx_map['symbol'].get(set_id.split('_')[0])
    elif set_id.startswith('1_') \
         or set_id.startswith('11_') \
         or set_id.startswith('12_'):
        return ras_chroms.get(set_id.split('_')[0])
    elif '_' in set_id:
        return set_id.split('_')[0]
    else:
        # If this point has been reached, need a new strategy to infer genes
        print('Unsure how to parse set ID ' + set_id + '\n')
        import pdb; pdb.set_trace()


def update_res(subres, infile, tx_map, min_freq=0.01):
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

    # Map categories onto cancer types & genes in subres
    for cancer in cancers:
        hits = df[cancer + '_AF'] >= min_freq
        for idx, vals in df.loc[hits, 'set_id gene'.split()].iterrows():
            set_id, gstr = vals.values
            for gene in gstr.split(','):
                if gene not in subres[cancer].keys():
                    continue
                # Convert DEL/AMP to simple set ID to avoid slightly different variant IDs
                clean_set_id = ''
                for sub_id in set_id.split('|'):
                    if any([sub_id.endswith(x) for x in '_DEL _AMP _DUP'.split()]):
                        sub_id = '_'.join([gene, set_id.split('_')[-1]])
                    clean_set_id = '|'.join([clean_set_id, sub_id])
                subres[cancer][gene].add(clean_set_id)

    # Return updated subres
    return subres


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutations', action='append', 
                        help='frequencies corresponding to individual mutations')
    parser.add_argument('--codons', action='append', 
                        help='frequencies corresponding to recurrent codon mutations')
    parser.add_argument('--burden-sets', action='append', 
                        help='frequencies corresponding to burden sets')
    parser.add_argument('--comutations', action='append', 
                        help='frequencies corresponding to comutation pairs')
    parser.add_argument('--ras-nonras-comut', action='append', 
                        help='frequencies corresponding to RAS + non-RAS ' + \
                        'comutation pairs')
    parser.add_argument('-t', '--transcript-info', required=True, help='.tsv ' + \
                        'mapping ENST:ENSG:symbol:length')
    parser.add_argument('-m', '--min-freq', default=0.01, type=float, help='Minimum ' + \
                        'frequency for a category to be retained per cancer type.')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Build dict for collecting results
    res = {cat : {cncr : {gene : set() for gene in ras_genes} for cncr in cancers} \
           for cat in category_descriptions.keys()}

    # Load transcript map
    tx_map = load_tx_map(args.transcript_info)
    
    # Load individual mutations
    for infile in args.mutations:
        res['mutations'] = \
            update_res(res['mutations'], infile, tx_map, args.min_freq)

    # Load recurrently mutated codons
    for infile in args.codons:
        res['codons'] = \
            update_res(res['codons'], infile, tx_map, args.min_freq)

    # Load burden sets
    for infile in args.burden_sets:
        res['burden'] = \
            update_res(res['burden'], infile, tx_map, args.min_freq)

    # Load comutation pairs
    for infile in args.comutations:
        res['comutations'] = \
            update_res(res['comutations'], infile, tx_map, args.min_freq)

    # Load RAS + non-RAS comutation pairs
    for infile in args.ras_nonras_comut:
        res['ras_nonras_comut'] = \
            update_res(res['ras_nonras_comut'], infile, tx_map, args.min_freq)

    # Open connection to --outfile
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    header_vals = ['Somatic Criteria']
    for cancer in cancers:
        for gene in ras_genes:
            header_vals.append('_'.join([cancer, gene]))
        header_vals.append(cancer + '_Union')
    outfile.write('\t'.join(header_vals + ['Total']) + '\n')

    # Summarize as table
    for category, vals in res.items():
        outvals = [category_descriptions[category]]
        total = 0
        for cancer in cancers:
            sub_union = set()
            for gene in ras_genes:
                outvals.append(len(vals[cancer][gene]))
                sub_union.update(vals[cancer][gene])
            outvals.append(len(sub_union))
            total += len(sub_union)
        outfile.write('\t'.join([str(x) for x in outvals + [total]]) + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

