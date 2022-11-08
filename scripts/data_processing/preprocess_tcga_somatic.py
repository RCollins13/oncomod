#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Subset TCGA somatic data to relevant tumors
"""


import argparse
import numpy as np
import pandas as pd
import pybedtools as pbt
from harmonize_tcga_samples import _parse_donor
from sys import stdout


def load_mutation_data(mc3_tsv, donors, purity_tsv):
    """
    Load & curate mutation data for donors of interest
    """

    cols_to_keep = 'Chromosome Start_Position Reference_Allele Tumor_Seq_Allele2 Hugo_Symbol Variant_Classification Tumor_Sample_Barcode HGVSc HGVSp_Short Transcript_ID Exon_Number t_depth t_alt_count Gene SYMBOL EXON VARIANT_CLASS FILTER'
    mut_df = pd.read_csv(mc3_tsv, sep='\t', header=0, low_memory=False,
                         usecols=cols_to_keep.split())

    # Add missing columns
    vid_cols = 'Chromosome Start_Position Reference_Allele ' + \
               'Tumor_Seq_Allele2 Tumor_Sample_Barcode'
    mut_df['VARIANT_CALL_ID'] = \
        mut_df[vid_cols.split()].astype(str).apply(lambda x: '_'.join(x), axis=1)
    mut_df['DONOR_ID'] = mut_df.Tumor_Sample_Barcode.astype(str).apply(_parse_donor)
    mut_df['ALLELE_FRACTION'] = (mut_df.t_alt_count / mut_df.t_depth)
    mut_df['TEST_TYPE'] = 'EXOME'
    na_cols = 'PATHOLOGIST_PATHOGENICITY MAX_GNOMAD_FREQUENCY PANEL_VERSION ' + \
              'CANONICAL_ENSEMBL_TSCP_ID CANONICAL_REF_SEQ_TSCP_ID'
    for col in na_cols.split():
        mut_df[col] = pd.NA
    mut_df['CANONICAL_EXON'] = mut_df.EXON.astype(str).apply(lambda x: x.split('/')[0])

    # Subset to PASS variants and donors in donor list
    mut_df = mut_df.loc[mut_df.DONOR_ID.isin(donors) & (mut_df.FILTER == 'PASS'), :]

    # Append purity data
    purity_df = pd.read_csv(purity_tsv, sep='\t', header=0)
    purity_df['Sample ID'] = purity_df['Sample ID'].apply(_parse_donor)
    purity_map = purity_df.set_index('Sample ID', drop=True).to_dict()['CPE']
    mut_df['TUMOR_PURITY'] = mut_df.DONOR_ID.map(purity_map)

    # Remap column names
    mut_df.rename(columns={'Chromosome' : 'CHROMOSOME',
                           'Start_Position' : 'POSITION',
                           'Reference_Allele' : 'REF_ALLELE',
                           'Tumor_Seq_Allele2' : 'ALT_ALLELE',
                           'SYMBOL' : 'CANONICAL_GENE',
                           'Hugo_Symbol' : 'HARMONIZED_HUGO_GENE_NAME',
                           'HGVSp_Short' : 'HARMONIZED_PROTEIN_CHANGE',
                           'HGVSc' : 'HARMONIZED_CDNA_CHANGE',
                           'Variant_Classification' : 'HARMONIZED_VARIANT_CLASS',
                           'Transcript_ID' : 'HARMONIZED_TRANSCRIPT_ID',
                           'Tumor_Sample_Barcode' : 'TUMOR_ID',
                           't_depth' : 'COVERAGE'},
                  inplace=True)

    # Add duplicate columns to match PROFILE data
    dup_col_suffixes = 'PROTEIN_CHANGE CDNA_CHANGE VARIANT_CLASS'
    for suf in dup_col_suffixes.split():
        mut_df['CANONICAL_' + suf] = mut_df['HARMONIZED_' + suf]

    # Return data frame sorted by chrom/pos and with reordered columns
    col_order = 'CHROMOSOME POSITION REF_ALLELE ALT_ALLELE VARIANT_CALL_ID DONOR_ID ' + \
                'TUMOR_ID CANONICAL_GENE HARMONIZED_HUGO_GENE_NAME ' + \
                'CANONICAL_PROTEIN_CHANGE HARMONIZED_PROTEIN_CHANGE ' + \
                'CANONICAL_CDNA_CHANGE HARMONIZED_CDNA_CHANGE ' + \
                'CANONICAL_VARIANT_CLASS HARMONIZED_VARIANT_CLASS CANONICAL_EXON ' + \
                'CANONICAL_ENSEMBL_TSCP_ID CANONICAL_REF_SEQ_TSCP_ID ' + \
                'HARMONIZED_TRANSCRIPT_ID PATHOLOGIST_PATHOGENICITY ALLELE_FRACTION ' + \
                'COVERAGE MAX_GNOMAD_FREQUENCY TUMOR_PURITY TEST_TYPE PANEL_VERSION'
    
    return mut_df.sort_values('CHROMOSOME POSITION'.split())[col_order.split()]


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mc3-tsv', help='Somatic mutation .tsv', required=True)
    parser.add_argument('--donors-list', help='List of donors to evaluate', required=True)
    parser.add_argument('--purity-tsv', help='Tumor purity .tsv', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to somatic ' +
                        'variation .tsv [default: stdout]')
    args = parser.parse_args()

    # Load donors of interest
    with open(args.donors_list) as fin:
        donors = [l.rstrip() for l in fin.readlines()]

    # Load mutation data and subset to samples of interest
    mut_df = load_mutation_data(args.mc3_tsv, donors, args.purity_tsv)
    
    # Write somatic mutation data to output file
    if args.outfile in '- stdout'.split():
        fout = stdout
    else:
        fout = args.outfile
    mut_df = mut_df.rename(columns={'CHROMOSOME' : '#CHROMOSOME'})
    mut_df.to_csv(fout, sep='\t', index=False, na_rep='.')


if __name__ == '__main__':
    main()

