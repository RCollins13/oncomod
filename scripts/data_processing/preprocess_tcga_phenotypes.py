#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Preprocess TCGA clinical data for relevant patients
"""


import argparse
import pandas as pd
import re
from harmonize_tcga_samples import _parse_donor


def load_id_map(tsv_in):
    """
    Load dataframe mapping various IDs
    """

    id_df = pd.read_csv(tsv_in, sep='\t')

    # Only map needed is array IDs to donor IDs
    return id_df.set_index('ARRAY_TYPED_ID')['#DONOR_ID'].to_dict()


def load_clinical(cdr_in, tcga_study_table_in, bmi_in, id_map):
    """
    Load & clean TCGA clinical data resource
    """

    # Load data
    keep_cols = 'bcr_patient_barcode type age_at_initial_pathologic_diagnosis gender ' + \
                'ajcc_pathologic_tumor_stage vital_status OS.time cause_of_death'
    main_df = pd.read_csv(cdr_in, usecols=keep_cols.split())

    # Subset to samples present in ID map
    main_df = main_df[main_df.bcr_patient_barcode.isin(id_map.values())]

    # Clean data
    main_df['CANCER_TYPE'] = main_df.type.map({'PAAD' : 'PDAC', 'COAD' : 'CRAD',
                                               'READ' : 'CRAD', 'SKCM' : 'SKCM'})
    main_df.ajcc_pathologic_tumor_stage = main_df.ajcc_pathologic_tumor_stage.str.replace('Stage ', '')
    for k, v in {'III' : '3', 'II' : '2', 'IV' : '4', 'I' : '1'}.items():
        main_df.ajcc_pathologic_tumor_stage = main_df.ajcc_pathologic_tumor_stage.str.replace(k, v)
    def _simplify_stage(sv):
        parseable = bool(re.search('^[0-4]', sv))
        if parseable:
            return sv[0]
        else:
            return pd.NA
    main_df['APPROX_STAGE'] \
        = main_df.ajcc_pathologic_tumor_stage.astype(str).apply(_simplify_stage)
    main_df['IS_ALIVE'] = main_df.vital_status.map({'Alive' : 1, 'Dead' : 0})
    main_df.loc[main_df.cause_of_death == '[Not Available]', 'cause_of_death'] = pd.NA

    # Add missing columns
    study_map = pd.read_csv(tcga_study_table_in, sep='\t', header=None).\
                   set_index(0, drop=True).to_dict()[1]
    main_df['PRIMARY_CANCER_DIAGNOSIS'] = main_df.type.map(study_map)
    bmi_df = pd.read_csv(bmi_in, sep='\t', header=None)
    bmi_df[0] = bmi_df[0].apply(_parse_donor)
    bmi_map = bmi_df.set_index(0, drop=True).to_dict()[1]
    main_df['BMI'] = main_df.bcr_patient_barcode.map(bmi_map)

    return main_df


def add_ancestry_info(main_df, ancestry_in, pcs_in, id_map):
    """
    Add ancestry & genetic PCs for all samples
    """

    # Add ancestry
    anc_map = pd.read_csv(ancestry_in, sep='\t').set_index('SampleID').\
                 Ancestry_assignment.to_dict()
    main_df['POPULATION'] = main_df.bcr_patient_barcode.map(anc_map)

    # Add genetic PCs
    pc_colnames = ['STUDY', 'SAMPLE'] + ['PC{}'.format(i+1) for i in range(10)]
    pc_df = pd.read_csv(pcs_in, sep=' ', header=None, names=pc_colnames).\
               drop('STUDY', axis=1)
    pc_df.SAMPLE = pc_df.SAMPLE.apply(_parse_donor)
    main_df = main_df.merge(pc_df, left_on='bcr_patient_barcode', 
                            right_on='SAMPLE', how='left')

    return main_df


def clean_output_df(main_df):
    """
    Clean up main patient data frame before writing to output file
    """

    new_col_names = {'bcr_patient_barcode' : 'DONOR_ID',
                     'age_at_initial_pathologic_diagnosis' : 'AGE_AT_DIAGNOSIS',
                     'OS.time' : 'DAYS_SURVIVED',
                     'cancer_type' : 'CANCER_TYPE',
                     'Predicted_Pop' : 'POPULATION',
                     'cause_of_death' : 'CAUSE_OF_DEATH',
                     'gender' : 'SEX',
                     'ajcc_pathologic_tumor_stage' : 'AJCC_STAGE'}
    missing_columns = 'DIAGNOSIS_DATE DIAGNOSIS_YEAR LAST_ALIVE_DATE ' + \
                      'BIOPSY_SITE BIOPSY_SITE_TYPE'
    drop_cols = 'type vital_status SAMPLE'

    main_df.rename(columns=new_col_names, inplace=True)
    main_df.drop(drop_cols.split(), axis=1, inplace=True)
    for col in missing_columns.split():
        main_df[col] = pd.NA

    col_order = 'DONOR_ID SEX BMI POPULATION CANCER_TYPE PRIMARY_CANCER_DIAGNOSIS ' + \
                'AJCC_STAGE APPROX_STAGE DIAGNOSIS_DATE AGE_AT_DIAGNOSIS ' + \
                'DIAGNOSIS_YEAR IS_ALIVE LAST_ALIVE_DATE DAYS_SURVIVED ' + \
                'CAUSE_OF_DEATH PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 ' + \
                'BIOPSY_SITE BIOPSY_SITE_TYPE'
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
    parser.add_argument('--cdr-csv', help='TCGA CDR .csv', required=True)
    parser.add_argument('--tcga-study-table', help='TCGA study code table', required=True)
    parser.add_argument('--bmi-tsv', help='TCGA BMI .tsv', required=True)
    parser.add_argument('--ancestry-tsv', help='Ancestry .tsv', required=True)
    parser.add_argument('--pcs-txt', help='Principal components .txt', required=True)
    parser.add_argument('--out-prefix', help='Prefix for output files', required=True)
    args = parser.parse_args()

    # Load ID map
    id_map = load_id_map(args.id_map_tsv)

    # Load & clean CDR and subset to samples of interest
    main_df = load_clinical(args.cdr_csv, args.tcga_study_table, 
                            args.bmi_tsv, id_map)

    # Add ancestry information
    main_df = add_ancestry_info(main_df, args.ancestry_tsv, args.pcs_txt, id_map)

    # Clean up output dataframe
    main_df = clean_output_df(main_df)

    # Write out full table of patient metadata for all cancer types
    main_df.rename(columns = {'DONOR_ID' : '#DONOR_ID'}, inplace=True)
    main_df.to_csv(args.out_prefix + 'ALL.sample_metadata.tsv.gz', 
                   index=False, sep='\t', na_rep='.')


if __name__ == '__main__':
    main()
