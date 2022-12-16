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
import pandas as pd
from sys import stdout, path
# TODO: uncomment the below before full commit
# path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
# from general_utils import load_tx_map


cancers = 'PDAC CRAD SKCM LUAD'.split()
ras_genes = 'NRAS HRAS KRAS'.split()
category_descriptions = \
    {'mutations' : '1. Frequent RAS mutations',
     'codons'    : '2. Recurrently mutated RAS codons',
     'burden'    : '3. Collapsed mutation sets'}


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


def update_res(subres, infile, tx_map):
    """
    Update category ID sets based on data in infile
    """

    # Load data from infile
    df = pd.read_csv(infile, sep='\t')

    # Attempt to infer gene for each category
    # TODO: implement this
    import pdb; pdb.set_trace()

    # Map categories onto cancer types & genes in subres
    # TODO: implement this

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
    parser.add_argument('--codons', action='append', nargs='+', 
                        help='frequencies corresponding to recurrent codon mutations')
    parser.add_argument('--burden-sets', action='append', nargs='+', 
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
        res['mutations'] = update_res(res['mutations'], infile, tx_map)





if __name__ == '__main__':
    main()

