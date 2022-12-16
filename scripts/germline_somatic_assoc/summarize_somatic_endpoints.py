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
# TODO: uncomment the below before full commit
# path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
# from general_utils import load_tx_map


# Define various variables used throughout the below functions
cancers = 'PDAC CRAD SKCM LUAD'.split()
ras_genes = 'NRAS HRAS KRAS'.split()
ras_chroms = {'1' : 'NRAS', '11' : 'HRAS', '12' : 'KRAS'}
category_descriptions = \
    {'mutations' : '1. Frequent RAS alterations',
     'codons'    : '2. Recurrently mutated RAS codons',
     'burden'    : '3. Collapsed mutation sets'}
tissue_map = {'PDAC' : 'pancreas',
              'CRAD' : 'colon',
              'SKCM' : 'skin',
              'LUAD' : 'lung'}


# TODO: remove the below function before full commit
def load_tx_map(tx_tsv):
    """
    Load --transcript-info as a dict
    """
    
    # Load data
    tx_df = pd.read_csv(tx_tsv, sep='\t', header=None)
    tx_df.columns = 'ENST ENSG symbol tx_len'.split()

    # Remove Ensembl version info from ENSG and ENST IDs
    for ecol in 'ENST ENSG'.split():
        tx_df[ecol] = tx_df[ecol].str.replace('.[0-9]+$', '', regex=True)

    # Map everything vs. ENST IDs
    tx_df.set_index('ENST', drop=True, inplace=True)
    
    return {'ENSG' : tx_df.ENSG.to_dict(),
            'ENSG_to_symbol' : tx_df.symbol.set_axis(tx_df.ENSG).to_dict(),
            'symbol' : tx_df.symbol.to_dict(),
            'tx_len' : tx_df.tx_len.to_dict()}


def update_res(subres, infile, tx_map, min_freq=0.01):
    """
    Update category ID sets based on data in infile
    """

    # Load data from infile
    df = pd.read_csv(infile, sep='\t')

    # Attempt to infer gene for each category
    if df.set_id[0].startswith('ENST'):
        df['gene'] = df.set_id.apply(lambda x: x.split('_')[0]).map(tx_map['symbol'])
    elif df.set_id[0].startswith('1_') \
         or df.set_id[0].startswith('11_') \
         or df.set_id[0].startswith('12_'):
        df['gene'] = df.set_id.apply(lambda x: x.split('_')[0]).map(ras_chroms)
    else:
        import pdb; pdb.set_trace()

    # Map categories onto cancer types & genes in subres
    for cancer in cancers:
        hits = df[cancer + '_AF'] > min_freq
        for idx, vals in df.loc[hits, 'set_id gene'.split()].iterrows():
            set_id, gene = vals.values
            # Convert DEL/AMP to simple set ID to avoid slightly different variant IDs
            if any([set_id.endswith(x) for x in '_DEL _AMP _DUP'.split()]):
                set_id = '_'.join([gene, set_id.split('_')[-1]])
            subres[cancer][gene].add(set_id)

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

