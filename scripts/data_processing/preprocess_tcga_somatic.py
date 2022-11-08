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


def load_mutation_data(mc3_tsv, donors, purity_tsv, key_genes=None):
    """
    Load & curate mutation data for donors of interest
    """

    cols_to_keep = 'Chromosome Start_Position Reference_Allele Tumor_Seq_Allele2 ' + \
                   'Hugo_Symbol Variant_Classification Tumor_Sample_Barcode ' + \
                   'HGVSc HGVSp_Short Transcript_ID Exon_Number t_depth ' + \
                   't_alt_count Gene SYMBOL EXON VARIANT_CLASS FILTER'
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

    # Subset to genes of interest, if optioned
    if key_genes is not None:
        mut_df = mut_df.loc[mut_df.Hugo_Symbol.isin(key_genes) | \
                            mut_df.SYMBOL.isin(key_genes)]

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
    
    return mut_df.sort_values('CHROMOSOME POSITION'.split())[col_order.split()], purity_map


def add_cna_data(mut_df, cna_bed, genes_gtf, donors, purity_map, key_genes=None):
    """
    Load CNA data, transform to gene-level results, and merge with mut_df
    """

    # Load CNA data
    cna_df = pd.read_csv(cna_bed, header=0, sep='\t').rename(columns={'sample' : 'SAMPLE'})
    cna_df['DONOR'] = cna_df.SAMPLE.apply(_parse_donor)
    cna_bt = pbt.BedTool.from_dataframe(cna_df.loc[cna_df.DONOR.isin(donors)])

    # Preprocess GTF
    def _preprocess_gtf(feature):
        feature.chrom = str(feature.chrom).replace('chr', '')
        pop_keys = [k for k in feature.attrs.keys() if k != 'gene_name']
        for key in pop_keys:
            feature.attrs.pop(key)
        return feature
    gene_bt = pbt.BedTool(genes_gtf).\
                  filter(lambda x: x.fields[2] == 'gene').\
                  each(_preprocess_gtf)
    if key_genes is not None:
        gene_bt = gene_bt.filter(lambda x: x.attrs.get('gene_name', None) in key_genes)

    # Intersect CNA and gene data and melt to donor-gene pairs compatible with mut_df
    tmp_df_names = 'CHROMOSOME cna_start cna_end cnv TUMOR_ID DONOR_ID gene_chrom ' + \
                   'source feature gene_start gene_end score strand skip attrs'
    hit_df = cna_bt.intersect(gene_bt, wb=True).to_dataframe(names=tmp_df_names.split())
    keep_cols = 'CHROMOSOME cnv TUMOR_ID DONOR_ID gene_start gene_end attrs'
    drop_cols = [k for k in tmp_df_names.split() if k not in keep_cols.split()]
    hit_df.drop(drop_cols, axis=1, inplace=True)
    hit_df['POSITION'] = \
        hit_df['gene_start gene_end'.split()].apply(np.nanmean, axis=1).astype(int)
    hit_df['REF_ALLELE'] = 'N'
    hit_df['ALT_ALLELE'] = hit_df.cnv.apply(lambda x: '<{}>'.format(x))
    for gcol in 'CANONICAL_GENE HARMONIZED_HUGO_GENE_NAME'.split():
        hit_df[gcol] = hit_df['attrs'].apply(lambda x: x.split('"')[1])
    hit_df['VARIANT_CALL_ID'] = \
        hit_df['DONOR_ID cnv CANONICAL_GENE'.split()].apply(lambda x: '_'.join(x), axis=1)
    for prefix in 'CANONICAL HARMONIZED'.split():
        for suffix in 'PROTEIN_CHANGE CDNA_CHANGE VARIANT_CLASS'.split():
            hit_df['{}_{}'.format(prefix, suffix)] = hit_df.cnv
    na_cols = 'CANONICAL_EXON CANONICAL_REF_SEQ_TSCP_ID CANONICAL_ENSEMBL_TSCP_ID ' + \
              'HARMONIZED_TRANSCRIPT_ID PATHOLOGIST_PATHOGENICITY ALLELE_FRACTION ' + \
              'COVERAGE MAX_GNOMAD_FREQUENCY PANEL_VERSION'
    for colname in na_cols.split():
        hit_df[colname] = pd.NA
    hit_df['TUMOR_PURITY'] = hit_df.DONOR_ID.map(purity_map)
    hit_df['TEST_TYPE'] = 'AFFY_SNP_6.0'
    
    # Add CNA data to mut_df, sort by position, and return
    out_df = pd.concat([hit_df.loc[:, mut_df.columns.tolist()], mut_df], ignore_index=True)
    return out_df.sort_values('CHROMOSOME POSITION'.split())


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mc3-tsv', help='Somatic mutation .tsv', required=True)
    parser.add_argument('--cna-bed', help='Copy number alteration .bed', required=True)
    parser.add_argument('--donors-list', help='List of donors to evaluate', required=True)
    parser.add_argument('--purity-tsv', help='Tumor purity .tsv', required=True)
    parser.add_argument('--genes-gtf', help='.gtf of gene annotations', required=True)
    parser.add_argument('--priority-genes', help='If provided, will subset all data ' + 
                        'to only variants in these genes [default: keep all genes]')
    parser.add_argument('-o', '--outfile', default='stdout', help='path to somatic ' +
                        'variation .tsv [default: stdout]')
    args = parser.parse_args()

    # Load donors of interest
    with open(args.donors_list) as fin:
        donors = list(set([l.rstrip() for l in fin.readlines()]))

    # Load genes of interest, if optioned
    if args.priority_genes is not None:
        with open(args.priority_genes) as fin:
            key_genes = list(set([l.rstrip() for l in fin.readlines()]))
    else:
        key_genes = None

    # Load mutation data and subset to samples of interest
    mut_df, purity_map = load_mutation_data(args.mc3_tsv, donors, 
                                            args.purity_tsv, key_genes)
    
    # Add CNA data to mut_df
    mut_df = add_cna_data(mut_df, args.cna_bed, args.genes_gtf, donors, 
                          purity_map, key_genes)
    
    # Write somatic mutation data to output file
    if args.outfile in '- stdout'.split():
        fout = stdout
    else:
        fout = args.outfile
    mut_df = mut_df.rename(columns={'CHROMOSOME' : '#CHROMOSOME'})
    mut_df.to_csv(fout, sep='\t', index=False, na_rep='.')


if __name__ == '__main__':
    main()

