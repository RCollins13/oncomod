#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Merge somatic variant data for plotting summary view of each RAS gene
"""


import argparse
import numpy as np
import os
import pandas as pd
from sys import path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from general_utils import load_tx_map


ras_genes = 'NRAS HRAS KRAS'.split()
cancers = 'PDAC CRAD LUAD SKCM'.split()


def load_coords(coords_in):
    """
    Load map of variant coordinates and split by cohort
    """

    df = pd.read_csv(coords_in, sep='\t')

    coords = {}
    for cohort, subdf in df.groupby('cohort'):
        sub_map = subdf.drop('cohort', axis=1).\
                        set_index('vid', drop=True).\
                        to_dict(orient='index')
        coords[cohort] = sub_map

    return coords


def load_set_map(set_map_in):
    """
    Load map of variant set memberships and split by cohort
    """

    df = pd.read_csv(set_map_in, sep='\t')
    df['vids'] = df.vids.str.split(',').map(lambda x: set(x))

    var_sets = {}
    for cohort, subdf in df.groupby('cohort'):
        sub_map = subdf.drop('cohort', axis=1).\
                        set_index('set_id', drop=True).\
                        vids.to_dict()
        var_sets[cohort] = sub_map

    return var_sets


def update_res_df(res_df, cohort, set_id, freqs, coords, var_sets, auto_add=False):
    """
    Update res_df for a single variant set in a single cohort 
    """

    # Get basic info
    vids = var_sets[cohort].get(set_id, set())

    # Attempt to find matching variant unless auto_add is specified
    if auto_add:
        
    else:
        # TODO: implement this
        import pdb; pdb.set_trace()

    return res_df



def unify_data(freqs, coords, var_sets, tx_map):
    """
    Collapse all relevant somatic mutation information into a single dense 
    table for plotting
    """

    cohorts = sorted(freqs.cohort.unique().tolist())

    # Link variants/sets between cohorts
    res_df = pd.DataFrame(columns='gene position csq'.split())
    for cohort, subdf in freqs.groupby('cohort'):

        # Add missing columns to res_df for cohort-specific info and cancer frequencies
        for suffix in 'set_id vids'.split():
            res_df[cohort + '_' + suffix] = [set()] * len(res_df)
        for cancer in cancers:
            for suffix in 'AN AC AF'.split():
                res_df['{}.{}_{}'.format(cohort, cancer, suffix)] = pd.NA

        # Update res_df for each variant set in serial
        for set_id in subdf.set_id.unique():
            res_df = update_res_df(res_df, cohort, set_id, freqs, coords, var_sets,
                                   auto_add=cohort == cohorts[0])

    return res_df


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--freqs', required=True, help='.tsv of variant frequencies, ' + \
                        'melted by cohort. Should follow output format from ' + \
                        'calc_mutation_freqs.py with appended left-most column ' + \
                        'indicating cohort.')
    parser.add_argument('--variant-coords', required=True, help='.tsv of coordinates ' + \
                        'for each variant, melted by cohort.')
    parser.add_argument('--variant-set-map', required=True, help='.tsv mapping ' + \
                        'collapsed variant sets to individual variant IDs.')
    parser.add_argument('-t', '--transcript-info', required=True, help='.tsv ' + \
                        'mapping ENST:ENSG:symbol:length')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Load frequencies and unify rows
    freqs = pd.read_csv(args.freqs, sep='\t')

    # Load map of variant coords
    coords = load_coords(args.variant_coords)
    
    # Load map of variant set memberships
    var_sets = load_set_map(args.variant_set_map)

    # Load transcript information
    tx_map = load_tx_map(args.transcript_info)

    # Unify all data
    res_df = unify_data(freqs, coords, var_sets, tx_map)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = args.outfile
    res_df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()

