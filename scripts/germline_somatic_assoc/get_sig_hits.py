#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Collect significant germline-somatic associations from one or more summary statistics files
"""


import argparse
import pandas as pd
from scipy.stats import norm
from sys import stdout


def load_ss(ss_tsv):
    """
    Load & process an input summary statistics .tsv
    """

    # Load data
    df = pd.read_csv(ss_tsv, sep='\t').rename(columns = {'#somatic' : 'somatic'})

    # Compute summary fields
    df['somatic_freq'] = df.somatic_AC / df.samples
    df['mutant_germline_CAF'] = df['yes_somatic.germline_AC'] / (2 * df.somatic_AC)
    df['control_germline_CAF'] = df['no_somatic.germline_AC'] / (2 * (df.samples - df.somatic_AC))

    # Drop unnecessary columns
    cols_to_drop = 'somatic_AC yes_somatic.germline_AC no_somatic.germline_AC ' + \
                   'z chisq model EUR_only'
    df.drop(columns = cols_to_drop.split(), inplace=True)

    return df


def get_sig_pairs(ss, p_cutoff=10e-8):
    """
    Extract all significant germline-somatic pairs from a dict of sumstats
    """

    sig_pairs = {}

    for df in ss.values():
        for vals in df[df.p <= p_cutoff].iterrows():
            sid = vals[1].somatic
            gid = vals[1].germline
            if sid not in sig_pairs.keys():
                sig_pairs[sid] = set([gid])
            else:
                sig_pairs[sid].add(gid)
    
    return sig_pairs


def format_output_res(ss, sig_pairs):
    """
    Format an output table of significant results
    """

    cohort_fields = 'samples somatic_freq mutant_caf control_caf beta p'.split()
    res_columns = 'somatic germline'.split()
    for cname in ss.keys():
        res_columns += ([cname + '_' + f for f in cohort_fields])

    res = pd.DataFrame(columns = res_columns)

    for sid, gids in sig_pairs.items():
        for gid in gids:
            out_vals = [sid, gid]
            for df in ss.values():
                in_vals = df[(df.somatic == sid) & (df.germline == gid)].squeeze()
                cohort_vals = [in_vals.samples, in_vals.somatic_freq, 
                               in_vals.mutant_germline_CAF,
                               in_vals.control_germline_CAF]
                ci_delta = norm.ppf(0.975) * in_vals.beta_SE
                lo_parts = [in_vals.beta, in_vals.beta - ci_delta, 
                            in_vals.beta + ci_delta]
                cohort_vals.append('{} [{}, {}]'.format(*[round(x, 2) for x in lo_parts]))
                cohort_vals.append(in_vals.p)
                out_vals += cohort_vals
            res.loc[len(res)] = out_vals

    # Sort by best P-value before returning
    p_cols = [cname + '_p' for cname in ss.keys()]
    res['best_p'] = res.loc[:, p_cols].min(axis=1)
    return res.sort_values('best_p').drop(columns='best_p')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--sumstats', action='append', 
                        help='path to one or more summary statistics .tsvs')
    parser.add_argument('-n', '--cohort-name', action='append', 
                        help='cohort names corresponding to --sumstats. Must be ' +
                        'supplied in the same order as --sumstats.')
    parser.add_argument('-p', '--p-cutoff', default=10e-8, type=float, 
                        help='P-value threshold to consider significant. ' + 
                        '[default: %default]')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Process each --sumstats input file
    ss = {cname : load_ss(ss_tsv) for ss_tsv, cname in zip(args.sumstats, args.cohort_name)}

    # Merge results for each significant hit
    sig_pairs = get_sig_pairs(ss, args.p_cutoff)
    res = format_output_res(ss, sig_pairs)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('\t\t' + ''.join([cname + '\t' * 6 for cname in ss.keys()]) + '\n')
    fields_per_cohort = ['Samples', 'Somatic Freq.', 'Mutant Germ. CAF', 
                         'Control Germ. CAF', 'Log Odds', 'P']
    outfile.write('Somatic\tGermline\t' + \
                  '\t'.join(fields_per_cohort * len(ss)) + '\n')
    res.to_csv(outfile, mode='a', header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()

