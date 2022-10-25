#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Preprocess GTEx v8 metasoft cross-tissue eQTL maps
Expected input: https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz
"""


import argparse
import pandas as pd
import numpy as np
import subprocess
from sys import stdout


def load_ensg_map(ensg_to_symbol_map):
    """
    Load map of ENSG-to-symbol conversions from input
    """

    with open(ensg_to_symbol_map, 'r') as tsvin:
        lines = [l.rstrip().split('\t') for l in tsvin.readlines()]
        return {ensg : gene for ensg, gene in lines}


def subset_to_sig_QTLs(raw_df):
    """
    Subset input QTL data to those with a Bonferroni-significant effect in at
    least one tissue
    """

    pval_df = raw_df[[k for k in raw_df.columns if k.startswith('pval_')]]
    best_p = pval_df.apply(np.nanmin, axis=1)
    best_p_is_gw = best_p <= (0.05 / raw_df['#STUDY'])

    return raw_df[best_p_is_gw]


def simplify_output(sig_df, ensg_map):
    """
    Simplify output data
    """

    # Get lists of gw & nominal sig tissues per variant
    pval_df = sig_df[[k for k in sig_df.columns if k.startswith('pval_')]]
    def __get_sig_colnames(data, pmin=0, pmax=1):
        sig = (data >= pmin) & (data < pmax)
        return [k.replace('pval_', '') for k in data.index[sig]]
    gw_sig_tissues = pval_df.apply(__get_sig_colnames, axis=1, pmax=10e-8)
    n_gw_sig = gw_sig_tissues.apply(len)
    nom_sig_tissues = pval_df.apply(__get_sig_colnames, axis=1, pmin=10e-8, pmax=0.05)
    n_nom_sig = nom_sig_tissues.apply(len)
    not_sig_tissues = pval_df.apply(__get_sig_colnames, axis=1, pmin=0.05)
    n_not_sig = not_sig_tissues.apply(len)

    # Fill empty cells with '.' for convenience
    def __fill_empties(data):
        if len(data) == 0:
            data.append('.')
        return data
    gw_sig_tissues = gw_sig_tissues.apply(__fill_empties)
    nom_sig_tissues = nom_sig_tissues.apply(__fill_empties)

    # Assign values to sig_df
    sig_df = sig_df.assign(n_gw_sig = n_gw_sig,
                           n_nom_sig = n_nom_sig,
                           n_not_sig = n_not_sig,
                           gw_sig_tissues = gw_sig_tissues.str.join(';'),
                           nom_sig_tissues = nom_sig_tissues.str.join(';'))

    # Extract variant information from RSID
    var_info = sig_df.RSID.str.split("_", expand=True)
    var_info.columns = '#chrom pos ref alt gene'.split()
    def __convert_gene(gstr, ensg_map):
        for ensg, gene in ensg_map.items():
            if ensg.split('.')[0] in gstr:
                return gene
    var_info = var_info.assign(gene = var_info.gene.apply(__convert_gene, ensg_map=ensg_map))

    # Reorder final columns for output df
    full_df = pd.concat([var_info, sig_df], axis=1)
    cols_to_keep = '#chrom pos ref alt gene n_gw_sig n_nom_sig n_not_sig ' + \
                   'gw_sig_tissues nom_sig_tissues'
    return full_df.loc[:, cols_to_keep.split()]


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metasoft_in', help='Input GTEx metasoft .tsv')
    parser.add_argument('--ensg-map', help='Two-column .tsv mapping ' +
                        'ENSG IDs to gene symbols', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to output ' +
                        'BED file [default: stdout]')
    parser.add_argument('--verbose-outfile', help='path to verbose output file ' +
                        '[default: no verbose outfile]')
    parser.add_argument('-z', '--bgzip', action='store_true', help='bgzip outfile')
    args = parser.parse_args()

    # Load input data as pd.DataFrame
    raw_df = pd.read_csv(args.metasoft_in, sep='\t')

    # Subset to QTLs with genome-wide significant effect in at least one tissue
    sig_df = subset_to_sig_QTLs(raw_df)
    if args.verbose_outfile is not None:
        sig_df.to_csv(args.verbose_outfile, sep='\t')

    # Load ENSG-to-symbol map as dict
    ensg_map = load_ensg_map(args.ensg_map)

    # Reformat simple output
    out_df = simplify_output(sig_df, ensg_map)

    # Write to outfile
    if args.outfile in '- stdout'.split():
        fout = stdout
        fout_path = 'stdout'
        stream_out = True
    else:
        stream_out = False
        if args.outfile.endswith('.gz'):
            gzip = True
            fout_path = args.outfile.replace('.gz', '')
        else:
            gzip = args.bgzip
            fout_path = args.outfile
        fout = open(fout_path, 'w')
    out_df.to_csv(fout, sep='\t', index=False)

    # Compress with bgzip, if optioned
    if gzip and not stream_out:
        subprocess.run(['bgzip', '-f', fout_path])


if __name__ == '__main__':
    main()
