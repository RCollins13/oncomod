#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Summarize number of germline collapsed burden sets to test for each cancer type & gene
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
tissue_map = {'PDAC' : 'pancreas',
              'CRAD' : 'colon',
              'SKCM' : 'skin',
              'LUAD' : 'lung'}


def update_res(res, infile, min_ac=10):
    """
    Update category ID sets based on data in infile
    """

    # Load data from infile
    df = pd.read_csv(infile, sep='\t')

    # Attempt to infer gene for each category
    df['gene'] = df.set_id.apply(lambda x: x.split('_')[0])

    # Check for tissue-specific annotations, and set frequencies for non-matched
    # cancers to zero to skip the counting step
    for cancer, tissue in tissue_map.items():
        tissue_rows = df.set_id.str.contains(tissue)
        if tissue_rows.any():
            other_ac_cols = [c for c in df.columns if c.endswith('_AC') and cancer not in c]
            df.loc[tissue_rows, other_ac_cols] = 0

    # Map categories onto cancer types & genes in res
    for cancer in cancers:
        hits = df[cancer + '_AC'] >= min_ac
        for idx, vals in df.loc[hits, 'set_id gene'.split()].iterrows():
            set_id, gstr = vals.values
            for gene in gstr.split(','):
                if gene not in res[cancer].keys():
                    continue
                res[cancer][gene].add(set_id)

    # Return updated res
    return res


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--burden-sets', action='append', 
                        help='frequencies corresponding to burden sets')
    parser.add_argument('-m', '--min-ac', default=10, type=int, help='Minimum ' + \
                        'AC for a category to be retained per cancer type.')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    parser.add_argument('-P', '--out-prefix', help='prefix for output lists of ' +
                        'germline sets to test [default: %default]',
                        type=str, default='./')
    args = parser.parse_args()

    # Build dict for collecting results
    res = {cncr : {gene : set() for gene in ras_genes} for cncr in cancers}
    
    # Load burden sets
    for infile in args.burden_sets:
        res = update_res(res, infile, args.min_ac)

    # Output lists of somatic endpoints per gene & cancer type
    for cancer in cancers:
        for gene in ras_genes:
            fout = open('{}{}.{}.germline_sets.tsv'.format(args.out_prefix, cancer, gene), 'w')
            for val in res[cancer][gene]:
                fout.write(val + '\n')
            fout.close()

    # Open connection to --outfile
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    header_vals = []
    for cancer in cancers:
        for gene in ras_genes:
            header_vals.append('_'.join([cancer, gene]))
        header_vals.append(cancer + '_Union')
    outfile.write('\t'.join(header_vals + ['Total']) + '\n')

    # Summarize as table
    total = 0
    outvals = []
    for cancer, vals in res.items():
        sub_union = set()
        for gene in ras_genes:
            outvals.append(len(vals[gene]))
            sub_union.update(vals[gene])
        outvals.append(len(sub_union))
        total += len(sub_union)
    outfile.write('\t'.join([str(x) for x in outvals + [total]]) + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

