#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Compute frequency of aribtrary set(s) of variant(s) from an allele dosage matrix
"""


import argparse
import numpy as np
import os
import pandas as pd
from sys import stdout, path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from freq_utils import load_variant_sets, load_cancer_map, load_dosage_matrix


def calc_af(var_sets, ad_df, cancer_map, max_an=2):
    """
    Compute AC, AN, and AF for each cancer type every entry in var_sets based on ad_df
    """

    # Get basic info for each cancer type
    cancers = sorted(list(set(cancer_map.values())))
    csamp_dict = {c : set() for c in cancers}
    for sample, cancer in cancer_map.items():
        csamp_dict[cancer].add(sample)

    # Prepare to collect results in dict keyed by variant set ID
    res_dict = {}

    # Compute values for each variant set in serial
    for sid, vids in var_sets.items():

        # First, check that all constituent variant IDs are present in ad_df
        if not all([v in ad_df.index for v in vids]):
            msg = 'not all variant IDs from {} present in --dosage-tsv'
            exit(msg.format(sid))

        # Second, check that variant set doesn't already exist in res_dict
        if sid in res_dict.keys():
            msg = 'variant set {} appears to be duplicated in --sets-tsv'
            exit(msg.format(sid))

        # If the above checks pass, then process each cancer type in serial
        res_dict[sid] = {}
        for cancer in cancers:
            vals = ad_df.loc[vids, ad_df.columns.isin(csamp_dict[cancer])].max()
            AN = int((~vals.isna()).sum() * max_an)
            AC = int(np.nansum(vals))
            if AN > 0:
                AF = AC / AN
            else:
                AF = np.nan
            res_dict[sid][cancer + '_AN'] = AN
            res_dict[sid][cancer + '_AC'] = AC
            res_dict[sid][cancer + '_AF'] = AF

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
    parser.add_argument('--sets-tsv', required=True, help='.tsv of variant set(s) ' +
                        'to evaluate. First and last columns must be set ID and ' + 
                        'comma-delimited list of variant IDs in set.')
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
    var_sets = load_variant_sets(args.sets_tsv)

    # Load dict mapping sample IDs to cancer type
    cancer_map = load_cancer_map(args.sample_metadata)

    # Load dosage matrix as pd.DataFrame
    ad_df = load_dosage_matrix(args.dosage_tsv, cancer_map.keys())

    # Compute AC, AN, and AF for every variant set by cancer type
    res_df = calc_af(var_sets, ad_df, cancer_map, args.max_an)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = args.outfile
    res_df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()

