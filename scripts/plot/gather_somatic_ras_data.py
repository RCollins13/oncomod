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
import re
from sys import path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from general_utils import load_tx_map


cancers = 'PDAC CRAD LUAD'.split()
aa_map = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", 
          "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", 
          "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", 
          "Tyr": "Y", "Val": "V"}


def load_coords(coords_in):
    """
    Load map of variant coordinates and split by cohort
    """

    df = pd.read_csv(coords_in, sep='\t')

    coords = {}
    for cohort, subdf in df.groupby('cohort'):
        try:
            sub_map = subdf.drop('cohort', axis=1).\
                            set_index('vid', drop=True).\
                            to_dict(orient='index')
        except:
            import pdb; pdb.set_trace()
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


def update_res_df(res_df, cohort, set_id, freqs, coords, var_sets, 
                  tx_map, auto_add=False):
    """
    Update res_df for a single variant set in a single cohort 
    """

    # Get basic info
    vids = var_sets[cohort].get(set_id, set())
    chrom = list(vids)[0].split('_')[0]
    all_pos = [int(v.split('_')[1]) for v in vids if 'AMP' not in v and 'DEL' not in v]
    if len(all_pos) > 0:
        pos = int(np.round(np.nanmean(all_pos), 0))
    else:
        pos = pd.NA
    if set_id.startswith('ENS'):
        sid_parts = set_id.split('_')
        enst = sid_parts[0]
        ensg = tx_map['ENSG'][enst]
        gene = tx_map['ENSG_to_symbol'][ensg]
        csq_wPrefix = sid_parts[1]
        if csq_wPrefix in 'AMP DEL'.split():
            alt = csq = csq_wPrefix
            codon = pd.NA
            ref = 'CN2'
        else:
            if len(sid_parts) > 2:
                csq = re.sub('^p\.', '', '_'.join(sid_parts[1:]))
            else:
                csq = re.sub('^p\.', '', csq_wPrefix)
            ref = re.split('[0-9]+', csq)[0]
            try:
                ref = aa_map[ref]
            except:
                pass
            alt = re.split('[0-9]+', csq, maxsplit=1)[-1]
            try:
                alt = aa_map[alt]
            except:
                pass
            codon = int(re.split('[A-z=]+', csq.split('_')[0])[1])

    else:
        gene = pd.NA
        codon = pd.NA
        alleles = set_id.split('_')[2:]
        if any([len(a) > 1 for a in alleles]):
            csq = 'other_InDel'
            ref, alt = alleles
        elif 'AMP' in set_id:
            ref = 'CN2'
            alt = csq = 'AMP'
        elif 'DEL' in set_id:
            ref = 'CN2'
            alt = csq = 'DEL'
        else:
            csq = 'other_SNV'
            ref, alt = alleles

    # Make small updates for CNAs
    if any([a in set_id for a in 'AMP DEL'.split()]):
        gene = chrom
        chrom = pd.NA

    # Attempt to find matching variant unless auto_add is specified
    freq_idx = (freqs.cohort == cohort) & (freqs.set_id == set_id)
    prior_idx = ((res_df.chrom == chrom) & (res_df.position == pos) & (res_df.alt == alt)) | \
                ((res_df.gene == gene) & (res_df.csq == csq))
    if auto_add or not prior_idx.any():
        new_vals = [gene, chrom, pos, csq, codon, ref, alt, set_id, ','.join(sorted(vids))]
        new_val_idxs = 'gene chrom position csq codon ref alt'.split()
        for suffix in 'set_id vids'.split():
            new_val_idxs.append(cohort + '_' + suffix)
        for suffix in freqs.columns[2:]:
            new_vals.append(freqs.loc[freq_idx, suffix].values[0])
            new_val_idxs.append(cohort + '_' + suffix)
        new_row = pd.Series(new_vals, index=new_val_idxs)
        res_df = res_df.append(new_row, ignore_index=True)
        
    else:
        res_df.loc[prior_idx, cohort + '_set_id'] = set_id
        res_df.loc[prior_idx, cohort + '_vids'] = ','.join(sorted(vids))
        for suffix in freqs.columns[2:]:
            newval = freqs.loc[freq_idx, suffix].values[0]
            res_df.loc[prior_idx, cohort + '_' + suffix] = newval

    return res_df


def unify_data(freqs, coords, var_sets, tx_map):
    """
    Collapse all relevant somatic mutation information into a single dense 
    table for plotting
    """

    cohorts = sorted(freqs.cohort.unique().tolist())

    # Link variants/sets between cohorts
    res_df = pd.DataFrame(columns='gene chrom position csq codon ref alt'.split())
    for cohort, subdf in freqs.groupby('cohort'):

        # Add missing columns to res_df for cohort-specific info and cancer frequencies
        for suffix in 'set_id vids'.split():
            res_df[cohort + '_' + suffix] = [pd.NA] * len(res_df)
        for cancer in cancers:
            for suffix in 'AN AC AF'.split():
                res_df['{}_{}_{}'.format(cohort, cancer, suffix)] = pd.NA

        # Update res_df for each variant set in serial
        for set_id in subdf.set_id.unique():
            res_df = update_res_df(res_df, cohort, set_id, freqs, coords, var_sets,
                                   tx_map, auto_add=cohort == cohorts[0])

    # Only retain variants with an annotated coding consequence
    res_df = res_df[~res_df.csq.str.contains('other_')]

    # Given that we are only retaining coding variants, it is ~safe to assume
    # that variants missing from each cohort were not observed (rather than)
    # artificially missing due to sequencing/technical factors
    # We will fill these NA columns to reflect AF~0
    for colname in res_df.columns[res_df.columns.str.endswith('_AN')]:
        res_df[colname] = res_df[colname].fillna(res_df[colname].max())
    for suffix in '_AC _AF'.split():
        for colname in res_df.columns[res_df.columns.str.endswith(suffix)]:
            res_df[colname] = res_df[colname].fillna(0)

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

    # Sort data by chrom, gene, codon, pos, alt
    res_df.sort_values('chrom gene codon position alt'.split(), inplace=True)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = args.outfile
    res_df.to_csv(outfile, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

