#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Merge LOHGIC predictions with OncoPanel mutation annotations
"""


import argparse
import pandas as pd


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--lohgic-tsv', help='LOHGIC output .tsv', 
                        required=True)
    parser.add_argument('--oncdrs-csv', help='OncDRS GENOMIC_MUTATION_RESULTS.csv',
                        required=True)
    parser.add_argument('-o', '--outfile', help='path to output file',
                        required=True)
    args = parser.parse_args()

    # Load LOHGIC results
    lohgic = pd.read_csv(args.lohgic_tsv, sep='\t')

    # Load mutations and retain relevant columns
    mut = pd.read_csv(args.oncdrs_csv, sep=',')
    drop_cols = 'DR_REQUEST_SEQ PATIENT_ID UNIQUE_SAMPLE_ID BEST_EFF_GENE BEST_EFF_GENE_ALIAS ' + \
                'BEST_EFF_PROTEIN_CHANGE BEST_EFF_CDNA_CHANGE BEST_EFF_VARIANT_CLASS ' + \
                'BEST_EFF_EXON BEST_EFFECT_ENSEMBL_TSCP_ID BEST_EFFECT_REF_SEQ_TSCP_ID ' + \
                'BEST_EFF_STRAND ALLELE_FRACTION COVERAGE'
    mut.drop(columns=drop_cols.split(), inplace=True)
    mut.rename(columns = {'TUMOR_PURITY' : 'PATHOLOGIST_PURITY',
                          'VARIANT_CALL_ID' : 'VARIANT_ID'}, 
               inplace=True)

    # Merge LOHGIC and OncDRS data
    mdf = lohgic.merge(mut, how='left', on='SAMPLE_ACCESSION_NBR VARIANT_ID'.split())

    # Reorder rows/columns and write to outfile
    col_order = 'CHROMOSOME POSITION REF_ALLELE ALT_ALLELE VARIANT_ID ' + \
                'SAMPLE_ACCESSION_NBR DFCI_MRN VARIANT_TYPE CANONICAL_GENE ' + \
                'CANONICAL_PROTEIN_CHANGE CANONICAL_CDNA_CHANGE ' + \
                'CANONICAL_VARIANT_CLASS CANONICAL_EXON CANONICAL_ENSEMBL_TSCP_ID ' + \
                'CANONICAL_REF_SEQ_TSCP_ID CANONICAL_STRAND HARMONIZED_HUGO_GENE_NAME ' + \
                'HARMONIZED_TRANSCRIPT_ID HARMONIZED_CDNA_CHANGE ' + \
                'HARMONIZED_PROTEIN_CHANGE HARMONIZED_VARIANT_CLASS ' + \
                'MAX_GNOMAD_FREQUENCY PLOIDY VAF DEPTH PURITY PURITY_CI ' + \
                'PATHOLOGIST_PURITY lohgic_somatic_weight lohgic_germline_weight ' + \
                'lohgic_model_1 lohgic_weight_1 lohgic_model_2 lohgic_weight_2 ' + \
                'PATHOLOGIST_TIER PATHOLOGIST_PATHOGENICITY RAPIDHEME_KB_PATHOGENICITY ' + \
                'PRIMARY_CANCER_DIAGNOSIS BIOPSY_SITE_TYPE BIOPSY_SITE ' + \
                'ORIGINAL_PATH_DIAGNOSIS TEST_TYPE PANEL_VERSION ' + \
                'CRDR_LAST_UPD_DTTM RECORD_CREATE_DT'
    mdf.loc[:, col_order.split()].\
        sort_values(by='CHROMOSOME POSITION ALT_ALLELE DFCI_MRN'.split()).\
        to_csv(args.outfile, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    main()

