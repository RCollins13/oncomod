#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Preprocess DFCI-Profile clinical data for relevant patients
"""


import argparse
import pandas as pd
from sys import stdout


def load_id_map(id_map_tsv, vcf_ids_in=None):
    """
    Load mappings of DFCI MRN to PBP IDs
    """

    id_df = pd.read_csv(id_map_tsv, sep='\t')['DFCI_MRN PBP'.split()]
    id_df.set_index('DFCI_MRN', inplace=True)
    
    # Subset to VCF IDs, if optioned
    if vcf_ids_in is not None:
        with open(vcf_ids_in) as fin:
            vcf_ids = set([i.rstrip() for i in fin.readlines()])
            id_df = id_df[id_df.PBP.isin(vcf_ids)]

    return id_df.to_dict()['PBP']


def add_dx_info(dx_csv):
    """
    # Load diagnostic table, remap IDs to PBPs, and restrict to samples with PBP IDs
    """

    cols_to_keep = 'PBP AGE_AT_DIAGNOSIS_NBR SITE_DESCR HISTOLOGY_DESCR'.split()

    dx_df = pd.read_csv(dx_csv, sep=',')
    dx_df = dx_df[dx_df.BEHAVIOR_DESCR == 'MALIGNANT, PRIMARY']
    dx_df['PBP'] = dx_df.DFCI_MRN.map(id_map)

    return dx_df.loc[~dx_df.PBP.isna(), cols_to_keep]


def add_ancestry_info(main_df, ancestry_csv, id_map):
    """
    Add ancestry information to main dataframe
    """

    cols_to_keep = 'Predicted_Pop ' + ' '.join(['PC' + str(i+1) for i in range(10)])
    cols_to_keep = cols_to_keep.split()

    ancestry_df = pd.read_csv(ancestry_tsv, sep=',')
    ancestry_df.set_index(ancestry_df.MRN.map(id_map), inplace=True)

    return main_df.merge(ancestry_df.loc[:, cols_to_keep], 
                         left_on='PBP', right_index=True)


def add_health_history(main_df, hx_csv, id_map):
    """
    Extract BMI for patients with available info
    """

    cols_to_read = 'DFCI_MRN D_START_DT HEALTH_HISTORY_TYPE RESULTS'.split()

    hx_df = pd.read_csv(hx_csv, sep=',', usecols=cols_to_read)
    hx_df = hx_df[hx_df.HEALTH_HISTORY_TYPE == 'BMI']
    hx_df['PBP'] = hx_df.DFCI_MRN.astype(int).map(id_map)
    hx_df = hx_df[~hx_df.PBP.isna()]
    
    bmi_map = {}
    for PBP in hx_df.PBP.unique():
        import pdb; pdb.set_trace()
        # hx_df.loc[hx_df.PBP == PBP]


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--id-map-tsv', help='.tsv mapping between various IDs', required=True)
    parser.add_argument('--dx-csv', help='Patient diagnosis .csv', required=True)
    parser.add_argument('--ancestry-csv', help='Ancestry .csv', required=True)
    parser.add_argument('--hx-csv', help='Health history .csv', required=True)
    parser.add_argument('--vcf-ids', dest='vcf_ids_in', help='Flat one-column ' +
                        'file with PBP IDs present in SNP VCF. If provided, ' +
                        'will subset all outputs to only these samples.')
    # parser.add_argument('-o', '--outfile', default='stdout', help='path to output ' +
    #                     'BED file [default: stdout]')
    # parser.add_argument('--verbose-outfile', help='path to verbose output file ' +
    #                     '[default: no verbose outfile]')
    # parser.add_argument('-z', '--bgzip', action='store_true', help='bgzip outfile')
    args = parser.parse_args()

    # Load ID map
    id_map = load_id_map(args.id_map_tsv, args.vcf_ids_in)

    # Add diagnostic information
    main_df = add_dx_info(args.dx_csv)

    # Add ancestry information
    main_df = add_ancestry_info(main_df, args.ancestry_csv, id_map)

    # Add relevant health history
    main_df = add_health_history(main_df, args.hx_csv, id_map)


    # panc_hits = dx_df.SITE_DESCR.str.contains('PANCR') \
    #             & dx_df.HISTOLOGY_DESCR.str.contains('ADENOCARCINOMA')
    # mel_hits = dx_df.SITE_DESCR.str.contains('SKIN') \
    #             & dx_df.HISTOLOGY_DESCR.str.contains('MELANOMA')
    # crc_hits = dx_df.HISTOLOGY_DESCR.str.contains('ADENOCARCINOMA') & ((dx_df.SITE_DESCR.str.contains('COLON')) | (dx_df.SITE_DESCR.str.contains('RECTUM')))


if __name__ == '__main__':
    main()
