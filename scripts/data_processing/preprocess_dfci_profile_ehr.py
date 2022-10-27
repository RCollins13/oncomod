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
import numpy as np
import pandas as pd
from collections import Counter
from sys import stdout


def _convert_DT_to_year(DT):
    """
    Helper function to convert a DT-style date into a four-digit year
    """

    try:
        short_year = int(DT.split('-')[2])
    except:
        return None

    if np.isnan(short_year):
        return np.nan
    elif short_year < 25:
        return int(2000 + short_year)
    else:
        return int(1900 + short_year)


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


def add_dx_info(dx_csv, id_map):
    """
    # Load diagnostic table, remap IDs to PBPs, and restrict to samples with PBP IDs
    """

    cols_to_keep = 'DFCI_MRN DIAGNOSIS_DT D_DIAGNOSIS_DT AGE_AT_DIAGNOSIS_NBR ' + \
                   'SITE_DESCR HISTOLOGY_DESCR BEHAVIOR_DESCR'

    dx_df = pd.read_csv(dx_csv, sep=',', usecols=cols_to_keep.split())
    dx_df = dx_df[dx_df.BEHAVIOR_DESCR == 'MALIGNANT, PRIMARY']
    dx_df['PBP'] = dx_df.DFCI_MRN.map(id_map)
    dx_df.drop('DFCI_MRN BEHAVIOR_DESCR'.split(), axis=1, inplace=True)
    dx_df['DIAGNOSIS_YEAR'] = dx_df.DIAGNOSIS_DT.astype(str).apply(_convert_DT_to_year)

    return dx_df.loc[~dx_df.PBP.isna()]


def subset_cancer_types(main_df):
    """
    Subset patient dataframe to cancer types of interest
    """

    # Find patient IDs matching relevant tumor types
    panc_hits = main_df.SITE_DESCR.str.contains('PANCR') \
                & main_df.HISTOLOGY_DESCR.str.contains('ADENOCARCINOMA')
    mel_hits = main_df.SITE_DESCR.str.contains('SKIN') \
               & main_df.HISTOLOGY_DESCR.str.contains('MELANOMA')
    crc_hits = main_df.HISTOLOGY_DESCR.str.contains('ADENOCARCINOMA') \
               & ((main_df.SITE_DESCR.str.contains('COLON')) | \
                  (main_df.SITE_DESCR.str.contains('RECTUM')))

    # Label patients based on their tumor types
    main_df['cancer_type'] = None
    main_df.loc[panc_hits, 'cancer_type'] = 'PDAC'
    main_df.loc[mel_hits, 'cancer_type'] = 'SKCM'
    main_df.loc[crc_hits, 'cancer_type'] = 'CRAD'

    # Drop patients with 2-3 relevant tumor types
    panc_ids = list(set(main_df.PBP[panc_hits]))
    mel_ids = list(set(main_df.PBP[mel_hits]))
    crc_ids = list(set(main_df.PBP[crc_hits]))
    dup_ids = [pbp for pbp, n in Counter(panc_ids + mel_ids + crc_ids).items() if n > 1]
    main_df[~main_df.PBP.isin(dup_ids)]

    # Subset to patients from cancer types of interest
    main_df = main_df[(panc_hits | mel_hits | crc_hits)]

    # For patients with multiple diagnostic entries, keep single row 
    # corresponding to earliest diagnosis
    main_df = main_df.sort_values('AGE_AT_DIAGNOSIS_NBR').\
                      drop_duplicates('PBP', keep='first')

    return main_df


def add_ancestry_info(main_df, ancestry_csv, id_map):
    """
    Add ancestry information to main dataframe
    """

    cols_to_keep = 'Predicted_Pop MRN ' + \
                   ' '.join(['PC' + str(i+1) for i in range(10)])

    ancestry_df = pd.read_csv(ancestry_csv, sep=',', usecols=cols_to_keep.split())
    ancestry_df.drop_duplicates('MRN', keep='first', inplace=True)
    ancestry_df.set_index(ancestry_df.MRN.map(id_map), inplace=True)
    ancestry_df.drop('MRN', axis=1, inplace=True)

    return main_df.merge(ancestry_df, left_on='PBP', right_index=True)


def add_health_history(main_df, hx_csv, id_map):
    """
    Extract BMI for patients with available info
    """

    cols_to_read = 'DFCI_MRN D_START_DT HEALTH_HISTORY_TYPE RESULTS'.split()
    phenos_to_keep = ['BMI']

    # Clean history dataframe
    hx_df = pd.read_csv(hx_csv, sep=',', usecols=cols_to_read)
    hx_df = hx_df[hx_df.HEALTH_HISTORY_TYPE.isin(phenos_to_keep)]
    hx_df['PBP'] = hx_df.DFCI_MRN.astype(int).map(id_map)
    hx_df = hx_df[~hx_df.PBP.isna()].sort_values('D_START_DT')
    hx_df.drop_duplicates('PBP', keep='first', inplace=True)
    
    # Extract BMI values
    bmi_map = hx_df[hx_df.HEALTH_HISTORY_TYPE == 'BMI'].set_index('PBP').RESULTS.to_dict()
    bmi_map = {k : float(v) for k, v in bmi_map.items()}
    main_df['BMI'] = main_df.PBP.map(bmi_map)

    return main_df


def add_survival_info(main_df, survival_csv, id_map):
    """
    Compute and add survival info for all patients
    """

    cols_to_keep = 'HYBRID_DEATH_IND HYBRID_DEATH_DT D_HYBRID_DEATH_DT ' + \
                   'NDI_CAUSE_OF_DEATH_RECODE DERIVED_LAST_ALIVE_DATE ' + \
                   'D_DERIVED_LAST_ALIVE_DATE GENDER_NM DFCI_MRN'

    surv_df = pd.read_csv(survival_csv, sep=',', low_memory=False, 
                          usecols=cols_to_keep.split())
    surv_df.drop_duplicates('DFCI_MRN', keep='first', inplace=True)
    surv_df.set_index(surv_df.DFCI_MRN.map(id_map), inplace=True)
    surv_df.drop('DFCI_MRN', axis=1, inplace=True)
    surv_df = surv_df[~surv_df.index.isna()]
    surv_df['IS_ALIVE'] = surv_df.HYBRID_DEATH_IND.map({'Y' : 0, 'N' : 1})

    main_df = main_df.merge(surv_df, left_on='PBP', right_index=True)
    main_df['DAYS_SURVIVED'] = main_df.D_DERIVED_LAST_ALIVE_DATE - main_df.D_DIAGNOSIS_DT

    return main_df


def clean_output_df(main_df):
    """
    Clean up main patient data frame before writing to output file
    """

    new_col_names = {'DIAGNOSIS_DT' : 'DIAGNOSIS_DATE',
                     'AGE_AT_DIAGNOSIS_NBR' : 'AGE_AT_DIAGNOSIS',
                     'cancer_type' : 'CANCER_TYPE',
                     'Predicted_Pop' : 'POPULATION',
                     'NDI_CAUSE_OF_DEATH_RECODE' : 'CAUSE_OF_DEATH',
                     'DERIVED_LAST_ALIVE_DATE' : 'LAST_ALIVE_DATE',
                     'GENDER_NM' : 'SEX'}
    drop_cols = 'D_DIAGNOSIS_DT SITE_DESCR HISTOLOGY_DESCR HYBRID_DEATH_IND ' + \
                'D_HYBRID_DEATH_DT D_DERIVED_LAST_ALIVE_DATE HYBRID_DEATH_DT'

    main_df.rename(columns=new_col_names, inplace=True)
    main_df.drop(drop_cols.split(), axis=1, inplace=True)

    col_order = 'PBP SEX BMI POPULATION CANCER_TYPE DIAGNOSIS_DATE ' + \
                'AGE_AT_DIAGNOSIS DIAGNOSIS_YEAR IS_ALIVE LAST_ALIVE_DATE ' + \
                'DAYS_SURVIVED CAUSE_OF_DEATH ' + \
                'PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10'
    sort_cols = 'CANCER_TYPE AGE_AT_DIAGNOSIS POPULATION'

    return main_df[col_order.split()].sort_values(sort_cols.split())


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
    parser.add_argument('--survival-csv', help='Survival info .csv', required=True)
    parser.add_argument('--out-prefix', help='Prefix for output files', required=True)
    parser.add_argument('--vcf-ids', dest='vcf_ids_in', help='Flat one-column ' +
                        'file with PBP IDs present in SNP VCF. If provided, ' +
                        'will subset all outputs to only these samples.')
    args = parser.parse_args()

    # Load ID map
    id_map = load_id_map(args.id_map_tsv, args.vcf_ids_in)

    # Add diagnostic information
    main_df = add_dx_info(args.dx_csv, id_map)

    # Subset to patients with cancer types of interest for study
    main_df = subset_cancer_types(main_df)

    # Add ancestry information
    main_df = add_ancestry_info(main_df, args.ancestry_csv, id_map)

    # Add relevant health history
    main_df = add_health_history(main_df, args.hx_csv, id_map)

    # Add survival information
    main_df = add_survival_info(main_df, args.survival_csv, id_map)

    # Write out PBP IDs per cancer type
    for cancer in 'PDAC SKCM CRAD'.split():
        with open(args.out_prefix + cancer + '.samples.list', 'w') as fout:
            ids = set(main_df.PBP[main_df.cancer_type == cancer].to_list())
            for sample in ids:
                fout.write(sample + '\n')
    all_cancers_out = args.out_prefix + 'ALL' + '.samples.list'
    main_df.loc[:, 'PBP cancer_type'.split()].\
            to_csv(all_cancers_out, sep='\t', index=False, header=False)

    # Write out full table of patient metadata for all cancer types
    main_df = clean_output_df(main_df)
    main_df.to_csv(args.out_prefix + 'ALL.sample_metadata.tsv.gz', 
                   index=False, sep='\t', na_rep='.')


if __name__ == '__main__':
    main()
