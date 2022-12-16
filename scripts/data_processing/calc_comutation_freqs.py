#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Compute comutation frequency of pairs of aribtrary set(s) of variant(s) from an allele dosage matrix
"""


import argparse
import numpy as np
import os
import pandas as pd
from itertools import combinations
from sys import path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from freq_utils import load_variant_sets, load_cancer_map, load_dosage_matrix


def calc_af(var_sets, pairs, ad_df, cancer_map, max_an=2):
    """
    Compute AC, AN, and AF for each cancer type for all pairs based on ad_df
    """

    # Get basic info for each cancer type
    cancers = sorted(list(set(cancer_map.values())))
    csamp_dict = {c : set() for c in cancers}
    for sample, cancer in cancer_map.items():
        csamp_dict[cancer].add(sample)

    # Prepare to collect results in dict keyed by comutation ID
    res_dict = {}

    # Compute values for each pair in serial
    for sid1, sid2 in pairs:

        # Get basic info
        comut_id = '|'.join(['COMUT', sid1, sid2])
        vids1 = var_sets[sid1]
        vids2 = var_sets[sid2]

        # First, check that all constituent variant IDs are present in ad_df
        if not all([v in ad_df.index for v in set(vids1).union(set(vids2))]):
            msg = 'not all variant IDs from {} present in --dosage-tsv'
            stop(msg.format(comut_id))

        # Second, check that comutation pair doesn't already exist in res_dict
        if comut_id in res_dict.keys():
            msg = 'comutation pair {} appears to be duplicated; are these specified ' + \
                  'more than once in --sets-tsv?'
            stop(msg.format(comut_id))

        # If the above checks pass, then process each cancer type in serial
        res_dict[comut_id] = {'sid1' : sid1, 'sid2' : sid2}
        for cancer in cancers:
            vals1 = ad_df.loc[vids1, ad_df.columns.isin(csamp_dict[cancer])].max()
            vals2 = ad_df.loc[vids2, ad_df.columns.isin(csamp_dict[cancer])].max()
            vals = (vals1 > 0) & (vals2 > 0)
            AN = int((~vals.isna()).sum() * max_an)
            AC = int(np.nansum(vals))
            if AN > 0:
                AF = AC / AN
            else:
                AF = np.nan
            res_dict[comut_id][cancer + '_AN'] = AN
            res_dict[comut_id][cancer + '_AC'] = AC
            res_dict[comut_id][cancer + '_AF'] = AF

    # Collapse all results into pd.DataFrame
    res_df = pd.DataFrame.from_dict(res_dict, orient='index')
    res_df.insert(loc=0, column='set_id', value=res_df.index)

    return res_df


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sets-tsv', required=True, action='append',
                        help='.tsv of variant set(s) to evaluate. First and ' + \
                        'last columns must be set ID and comma-delimited list ' + \
                        'of variant IDs in set. May be specified multiple times, ' + \
                        'in which case all sets from all files will be used.')
    parser.add_argument('--dosage-tsv', required=True, help='allele dosage matrix ' +
                        '.tsv. Rows = variants. First column = variant ID, ' + 
                        'remaining columns = samples.')
    parser.add_argument('--sample-metadata', required=True, help='sample metadata .tsv')
    parser.add_argument('--max-an', default=2, type=int, help='Max AN for all ' +
                        'sites [default: 2]')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Load map of variant sets
    var_sets = {}
    for fin in args.sets_tsv:    
        var_sets.update(load_variant_sets(fin))

    # Load dict mapping sample IDs to cancer type
    cancer_map = load_cancer_map(args.sample_metadata)

    # Load dosage matrix as pd.DataFrame
    ad_df = load_dosage_matrix(args.dosage_tsv, cancer_map.keys())

    # Enumerate all pairs of sets
    pairs = [x for x in combinations(var_sets.keys(), 2)]

    # Compute AC, AN, and AF for all pairs by cancer type
    res_df = calc_af(var_sets, pairs, ad_df, cancer_map, args.max_an)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = args.outfile
    res_df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()

