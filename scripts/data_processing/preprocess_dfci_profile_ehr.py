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
import re
from collections import Counter


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


def _status_report(main_df, stage=None):
    """
    Print a simple status report of sample counts per cancer type
    """

    preface = '\nSample count'
    if stage is not None:
        preface += ' after ' + stage
    print(preface + ':')

    for cancer, count in main_df.CANCER_TYPE.value_counts().iteritems():
        print('{}: {:,}'.format(cancer, int(count)))


def load_primary_data(genomic_csv, id_map_tsv, vcf_ids_in=None, 
                      priority_ids=None, allfit_tsv=None):
    """
    Load basic information for samples to be included in study
    Subset to cancer types of interest
    Map DFCI MRNs to PBP IDs
    Select most relevant PBP ID for MRNs with multiple PBPs
    """

    # Load genomic information
    genomic_cols_to_keep = 'DFCI_MRN SAMPLE_ACCESSION_NBR PRIMARY_CANCER_DIAGNOSIS ' + \
                           'BIOPSY_SITE BIOPSY_SITE_TYPE PANEL_VERSION CANCER_TYPE ' + \
                           'TUMOR_PURITY HOMOPOLYMER_PER_MEGABASE ' + \
                           'MISMATCH_REPAIR_STATUS'
    genomic_df = pd.read_csv(genomic_csv, sep=',', low_memory=False, 
                             usecols=genomic_cols_to_keep.split())
    genomic_df.rename(columns={'SAMPLE_ACCESSION_NBR' : 'BL_ID'}, inplace=True)
    genomic_df.loc[genomic_df.TUMOR_PURITY < 0, 'TUMOR_PURITY'] = pd.NA
    genomic_df.TUMOR_PURITY = genomic_df.TUMOR_PURITY / 100

    # Update purity from allfit_tsv if provided
    if allfit_tsv is not None:
        allfit_df = pd.read_csv(allfit_tsv, sep='\t').\
                       rename(columns={'PURITY' : 'TUMOR_PURITY',
                                       'SAMPLE_ACCESSION_NBR' : 'BL_ID'}).\
                       set_index('BL_ID', drop=True)
        allfit_dat = allfit_df.TUMOR_PURITY.to_dict()
        allfit_overlaps = genomic_df.BL_ID.isin(allfit_dat.keys())
        genomic_df.loc[allfit_overlaps, 'TUMOR_PURITY'] \
            = genomic_df.BL_ID[allfit_overlaps].map(allfit_dat)

    # Drop MSI tumors
    msi_idx = (genomic_df.HOMOPOLYMER_PER_MEGABASE > 2) | \
              (genomic_df.MISMATCH_REPAIR_STATUS == 'Deficient (MMR-D / MSI-H)')
    genomic_df = genomic_df.loc[~msi_idx, :]
    genomic_df.drop(columns='HOMOPOLYMER_PER_MEGABASE MISMATCH_REPAIR_STATUS'.split(), 
                    inplace=True)

    # Load all ID mappings
    id_cols_to_keep = 'DFCI_MRN SAMPLE_ACCESSION_NBR PBP'
    id_df = pd.read_csv(id_map_tsv, sep='\t', low_memory=False, 
                        usecols=id_cols_to_keep.split())
    id_df.rename(columns={'SAMPLE_ACCESSION_NBR' : 'BL_ID'}, inplace=True)

    # Subset genomic information to tumors of interest
    # According to Sasha, imputation quality differs among panels as follows:
    # panel 2 > panel 1 > panel 3.1 > panel 3
    # If the same patient was run on multiple panels, keep the one with highest priority
    # where priority is defined first by --priority-ids (if provided) 
    # and second by panel version. After that, arbitrarily choose just one 
    # of the tumor biopsies to keep per patient
    if priority_ids is not None:
        with open(priority_ids) as fin:
            pids_list = set([l.rstrip() for l in fin.readlines()])
        genomic_df.loc[:, 'LOW_PRIORITY'] = \
            ~genomic_df.BL_ID.map(id_df.set_index('BL_ID').to_dict()['PBP']).\
                              isin(pids_list)
    else:
        genomic_df.loc[:, 'LOW_PRIORITY'] = True
    panel_priority = [2.0, 1.0, 3.1, 3.0]
    genomic_df.PANEL_VERSION = pd.Categorical(genomic_df.PANEL_VERSION, panel_priority)
    biopsy_priority = 'PRIMARY,LOCAL RECURRENCE,DIAGNOSIS PENDING,' + \
                      'DIAGNOSIS UNAVAILABLE,UNSPECIFIED,ANY/OTHER,' + \
                      'METASTATIC RECURRENCE,NOT APPLICABLE'
    genomic_df.BIOPSY_SITE_TYPE = pd.Categorical(genomic_df.BIOPSY_SITE_TYPE.astype(str), 
                                                 biopsy_priority.split(','))
    genomic_df.sort_values('LOW_PRIORITY PANEL_VERSION BIOPSY_SITE_TYPE'.split(), 
                           inplace=True)
    panc_hits = (genomic_df.CANCER_TYPE == 'Pancreatic Cancer') \
                & (genomic_df.PRIMARY_CANCER_DIAGNOSIS.str.lower().str.contains('adenocarcinoma'))
    crc_hits = (genomic_df.CANCER_TYPE == 'Colorectal Cancer') \
               & (genomic_df.PRIMARY_CANCER_DIAGNOSIS.str.lower().str.contains('adenocarcinoma'))
    lung_hits = (genomic_df.CANCER_TYPE.str.contains('Lung Cancer')) \
                 & (genomic_df.PRIMARY_CANCER_DIAGNOSIS.str.lower().str.contains('adenocarcinoma'))
    genomic_df = genomic_df[(panc_hits | crc_hits | lung_hits)]
    genomic_df.drop_duplicates('DFCI_MRN', keep='first', inplace=True)
    genomic_df.drop(columns='LOW_PRIORITY', inplace=True)

    # Subset ID linker table to BL IDs present in genomic df
    id_df = id_df[id_df.BL_ID.isin(genomic_df.BL_ID)]

    # Subset to VCF IDs, if optioned
    if vcf_ids_in is not None:
        with open(vcf_ids_in) as fin:
            vcf_ids = set([i.rstrip() for i in fin.readlines()])
            id_df = id_df[id_df.PBP.isin(vcf_ids)]

    # At this point, there should be a 1:1:1 mapping between MRNs, BL IDs, and PBP IDs
    # Now we can build ID maps between the various formats
    id_map = {'MRN_to_BL' : id_df.set_index('DFCI_MRN').to_dict()['BL_ID'],
              'BL_to_MRN' : id_df.set_index('BL_ID').to_dict()['DFCI_MRN'],
              'MRN_to_PBP' : id_df.set_index('DFCI_MRN').to_dict()['PBP'],
              'PBP_to_MRN' : id_df.set_index('PBP').to_dict()['DFCI_MRN'],
              'BL_to_PBP' : id_df.set_index('BL_ID').to_dict()['PBP'],
              'PBP_to_BL' : id_df.set_index('PBP').to_dict()['BL_ID']}

    # Add PBP IDs to genomic data
    genomic_df['PBP'] = genomic_df.BL_ID.map(id_map['BL_to_PBP'])

    # Samples *must* have PBP IDs to be used for downstream analysis
    return genomic_df[~genomic_df.PBP.isna()]


def add_dx_info(main_df, dx_csv):
    """
    # Load diagnostic information for samples with available info
    """

    # Load dx information
    dx_cols_to_keep = 'DFCI_MRN DIAGNOSIS_DT D_DIAGNOSIS_DT AGE_AT_DIAGNOSIS_NBR ' + \
                      'SITE_DESCR HISTOLOGY_DESCR BEHAVIOR_DESCR BEST_AJCC_STAGE_CD'
    dx_df = pd.read_csv(dx_csv, sep=',', usecols=dx_cols_to_keep.split())

    # Subset diagnostic information to tumors/patients of interest
    dx_df = dx_df[dx_df.DFCI_MRN.isin(main_df.DFCI_MRN) \
                  & (dx_df.BEHAVIOR_DESCR == 'MALIGNANT, PRIMARY')]
    panc_hits = dx_df.SITE_DESCR.str.contains('PANCR')
    crc_hits = ((dx_df.SITE_DESCR.str.contains('COLON')) | \
                (dx_df.SITE_DESCR.str.contains('RECTUM')))
    lung_hits = dx_df.SITE_DESCR.str.contains('LUNG')
    dx_df = dx_df[panc_hits | crc_hits | lung_hits]
    
    # For the subset of patients with multiple diagnoses, take the dx info for
    # the earliest diagnosed tumor
    dx_df = dx_df.sort_values('D_DIAGNOSIS_DT').drop_duplicates('DFCI_MRN', keep='first')

    # Simplify stage info
    for k, v in {'III' : '3', 'II' : '2', 'IV' : '4', 'I' : '1'}.items():
        dx_df.BEST_AJCC_STAGE_CD = dx_df.BEST_AJCC_STAGE_CD.str.replace(k, v)
    def _simplify_stage(sv):
        parseable = bool(re.search('^[0-4]', sv))
        if parseable:
            return sv[0]
        else:
            return pd.NA
    dx_df['APPROX_STAGE'] = dx_df.BEST_AJCC_STAGE_CD.astype(str).apply(_simplify_stage)

    # Clean & merge dx info into main df
    dx_df.drop('BEHAVIOR_DESCR'.split(), axis=1, inplace=True)
    dx_df['DIAGNOSIS_YEAR'] = dx_df.DIAGNOSIS_DT.astype(str).apply(_convert_DT_to_year)
    main_df.merge(dx_df, left_on='DFCI_MRN', right_on='DFCI_MRN')

    return main_df.merge(dx_df, on='DFCI_MRN', how='left')


def add_ancestry_info(main_df, ancestry_csv):
    """
    Add ancestry information to main dataframe
    """

    cols_to_keep = 'Predicted_Pop X.FID ' + \
                   ' '.join(['PC' + str(i+1) for i in range(10)])

    ancestry_df = pd.read_csv(ancestry_csv, sep=',', usecols=cols_to_keep.split())
    ancestry_df.rename(columns={'X.FID' : 'PBP'}, inplace=True)

    return main_df.merge(ancestry_df, on='PBP', how='left')


# def add_health_history(main_df, hx_csv):
#     """
#     Extract BMI for patients with available info
#     """

#     cols_to_read = 'DFCI_MRN D_START_DT HEALTH_HISTORY_TYPE RESULTS'.split()
#     phenos_to_keep = ['BMI']

#     # Clean history dataframe
#     hx_df = pd.read_csv(hx_csv, sep=',', usecols=cols_to_read)
#     hx_df = hx_df[hx_df.HEALTH_HISTORY_TYPE.isin(phenos_to_keep)]
#     hx_df = hx_df.sort_values('D_START_DT').drop_duplicates('DFCI_MRN', keep='first')

#     return main_df


def add_survival_info(main_df, survival_csv):
    """
    Compute and add survival info for all patients
    """

    cols_to_keep = 'HYBRID_DEATH_IND HYBRID_DEATH_DT D_HYBRID_DEATH_DT ' + \
                   'NDI_CAUSE_OF_DEATH_RECODE DERIVED_LAST_ALIVE_DATE ' + \
                   'D_DERIVED_LAST_ALIVE_DATE GENDER_NM DFCI_MRN'

    surv_df = pd.read_csv(survival_csv, sep=',', low_memory=False, 
                          usecols=cols_to_keep.split())
    surv_df.drop_duplicates('DFCI_MRN', keep='first', inplace=True)
    surv_df['IS_ALIVE'] = surv_df.HYBRID_DEATH_IND.map({'Y' : 0, 'N' : 1})

    main_df = main_df.merge(surv_df, on='DFCI_MRN', how='left')
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
                     'GENDER_NM' : 'SEX',
                     'BEST_AJCC_STAGE_CD' : 'AJCC_STAGE'}
    drop_cols = 'D_DIAGNOSIS_DT SITE_DESCR HISTOLOGY_DESCR HYBRID_DEATH_IND ' + \
                'D_HYBRID_DEATH_DT D_DERIVED_LAST_ALIVE_DATE HYBRID_DEATH_DT ' + \
                'PANEL_VERSION'

    main_df.rename(columns=new_col_names, inplace=True)
    main_df.drop(drop_cols.split(), axis=1, inplace=True)

    # Simplify cancer types into codes
    cancer_codes = {'Pancreatic Cancer' : 'PDAC',
                    'Colorectal Cancer' : 'CRAD',
                    'Non-Small Cell Lung Cancer' : 'LUAD'}
    main_df.CANCER_TYPE = main_df.CANCER_TYPE.map(cancer_codes)

    col_order = 'PBP SEX POPULATION CANCER_TYPE PRIMARY_CANCER_DIAGNOSIS ' + \
                'AJCC_STAGE APPROX_STAGE DIAGNOSIS_DATE AGE_AT_DIAGNOSIS ' + \
                'DIAGNOSIS_YEAR IS_ALIVE LAST_ALIVE_DATE DAYS_SURVIVED ' + \
                'CAUSE_OF_DEATH PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 ' + \
                'BIOPSY_SITE BIOPSY_SITE_TYPE TUMOR_PURITY'
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
    parser.add_argument('--genomic-csv', help='Tumor genomic .csv', required=True)
    parser.add_argument('--dx-csv', help='Patient diagnosis .csv', required=True)
    parser.add_argument('--ancestry-csv', help='Ancestry .csv', required=True)
    # parser.add_argument('--hx-csv', help='Health history .csv')
    parser.add_argument('--survival-csv', help='Survival info .csv', required=True)
    parser.add_argument('--allfit-tsv', help='Refined purity estimates from All-FIT')
    parser.add_argument('--out-prefix', help='Prefix for output files', required=True)
    parser.add_argument('--vcf-ids', dest='vcf_ids_in', help='Flat one-column ' +
                        'file with PBP IDs present in SNP VCF. If provided, ' +
                        'will subset all outputs to only these samples.')
    parser.add_argument('--priority-ids', help='List of higher-priority PBP IDs. ' +
                        'These will always be retained over other samples from ' +
                        'the same patient wherever possible.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Report ' +
                        'sample counts per cancer type at each step ' + 
                        '[default: no reporting]')
    args = parser.parse_args()

    # Load primary data for all samples
    main_df = load_primary_data(args.genomic_csv, args.id_map_tsv, 
                                args.vcf_ids_in, args.priority_ids,
                                args.allfit_tsv)
    if args.verbose:
        _status_report(main_df, stage="loading primary data")

    # Add diagnostic information
    main_df = add_dx_info(main_df, args.dx_csv)
    if args.verbose:
        _status_report(main_df, stage="adding diagnostic information")

    # Add ancestry information
    main_df = add_ancestry_info(main_df, args.ancestry_csv)
    if args.verbose:
        _status_report(main_df, stage="annotating ancestry information")

    # # Add relevant health history
    # main_df = add_health_history(main_df, args.hx_csv)
    # if args.verbose:
    #     _status_report(main_df, stage="adding health history")

    # Add survival information
    main_df = add_survival_info(main_df, args.survival_csv)
    if args.verbose:
        _status_report(main_df, stage="adding survival data")

    # Clean up columns for final outputs
    main_df = clean_output_df(main_df)

    # Write out PBP IDs per cancer type
    for cancer in 'PDAC CRAD LUAD'.split():
        with open(args.out_prefix + cancer + '.samples.list', 'w') as fout:
            ids = set(main_df.PBP.astype(str)[main_df.CANCER_TYPE == cancer].to_list())
            for sample in ids:
                fout.write(sample + '\n')
    all_cancers_out = args.out_prefix + 'ALL' + '.samples.with_cancer_labels.tsv'
    main_df.loc[:, 'PBP CANCER_TYPE'.split()].\
            to_csv(all_cancers_out, sep='\t', index=False, header=False)
    all_samples_out = args.out_prefix + 'ALL' + '.samples.list'
    main_df.PBP.to_csv(all_samples_out, sep='\t', index=False, header=False)

    # Write out full table of patient metadata for all cancer types
    main_df.rename(columns = {'PBP' : '#PBP'}, inplace=True)
    main_df.to_csv(args.out_prefix + 'ALL.sample_metadata.tsv.gz', 
                   index=False, sep='\t', na_rep='.')


if __name__ == '__main__':
    main()
