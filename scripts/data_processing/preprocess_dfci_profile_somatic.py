#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Subset DFCI-Profile somatic data to relevant tumors
"""


import argparse
import numpy as np
import pandas as pd
import pybedtools as pbt
from sys import stdout


def load_id_map(samples_in, id_map_tsv):
    """
    Load map of various sample IDs and subset to samples in input list
    """

    with open(samples_in) as fin:
        samples = set([s.rstrip() for s in fin.readlines()])

    id_df = pd.read_csv(id_map_tsv, sep='\t')
    id_df = id_df[id_df.PBP.isin(samples)]

    return id_df.set_index('SAMPLE_ACCESSION_NBR', drop=True).to_dict()['PBP']


def load_mutation_data(mutation_csv, id_map):
    """
    Load & curate mutation data for samples of interest
    """

    cols_to_keep = 'CHROMOSOME POSITION REF_ALLELE ALT_ALLELE VARIANT_CALL_ID ' + \
                   'SAMPLE_ACCESSION_NBR CANONICAL_GENE HARMONIZED_HUGO_GENE_NAME ' + \
                   'CANONICAL_PROTEIN_CHANGE HARMONIZED_PROTEIN_CHANGE ' + \
                   'CANONICAL_CDNA_CHANGE HARMONIZED_CDNA_CHANGE ' + \
                   'CANONICAL_VARIANT_CLASS HARMONIZED_VARIANT_CLASS CANONICAL_EXON ' + \
                   'CANONICAL_ENSEMBL_TSCP_ID CANONICAL_REF_SEQ_TSCP_ID ' + \
                   'HARMONIZED_TRANSCRIPT_ID PATHOLOGIST_PATHOGENICITY ' + \
                   'ALLELE_FRACTION COVERAGE MAX_GNOMAD_FREQUENCY TUMOR_PURITY ' + \
                   'TEST_TYPE PANEL_VERSION'
    mut_df = pd.read_csv(mutation_csv, sep=',', usecols=cols_to_keep.split(), 
                         low_memory=False)
    
    # Remap IDs and drop to samples of interest
    mut_df['PBP'] = mut_df.SAMPLE_ACCESSION_NBR.map(id_map)
    mut_df = mut_df.loc[~mut_df.PBP.isna()]
    mut_df.drop('SAMPLE_ACCESSION_NBR', axis=1, inplace=True)

    # Sort columns and rows
    ordered_cols = cols_to_keep.split()
    i = ordered_cols.index('SAMPLE_ACCESSION_NBR')
    ordered_cols = ordered_cols[:i] + ['PBP'] + ordered_cols[i+1:]

    return mut_df.sort_values('CHROMOSOME POSITION'.split())[ordered_cols]


def add_cna_data(cna_csv, mut_df, id_map, genes_gtf=None):
    """
    Add CNA data on tumor samples of interest
    """

    cols_to_keep = 'GENE CNV_TYPE_CD COPY_COUNT CNV_CALL_ID SAMPLE_ACCESSION_NBR ' + \
                    'HARMONIZED_HUGO_GENE_NAME TUMOR_PURITY TEST_TYPE PANEL_VERSION'
    cna_df = pd.read_csv(cna_csv, sep=',', usecols=cols_to_keep.split())

    # Remap IDs and drop to samples of interest
    cna_df['PBP'] = cna_df.SAMPLE_ACCESSION_NBR.map(id_map)
    cna_df = cna_df.loc[~cna_df.PBP.isna()]
    cna_df.drop('SAMPLE_ACCESSION_NBR', axis=1, inplace=True)

    # If genes_gtf provided, build a map of chrom:pos for each gene
    gene_map = {}
    if genes_gtf is not None:
        for feature in pbt.BedTool(genes_gtf):
            if feature.fields[2] == 'gene':
                gene = feature.attrs['gene_name']
                gene_map[gene] = [feature.chrom.replace('chr', ''), np.nanmean([feature.start, feature.end])]

    # Convert alt alleles
    cna_df['REF_ALLELE'] = 'N'
    def __build_alt(cnv_data):
        cn = cnv_data.COPY_COUNT
        call = cnv_data.CNV_TYPE_CD
        if pd.isna(cn) or cn == 2.0:
            if call in '1DEL 2DEL'.split():
                return '<DEL>'
            elif call in 'LA HA'.split():
                return '<AMP>'
        else:
            if cn < 2:
                return '<DEL:CN{}>'.format(int(cn))
            elif cn > 2:
                return '<AMP:CN{}>'.format(int(cn))
    cna_df['ALT_ALLELE'] = cna_df['COPY_COUNT CNV_TYPE_CD'.split()].apply(__build_alt, axis=1)

    # Label CNV type for convenience
    def __label_cnv_type(alt):
        return alt.split(':')[0].replace('<', '').replace('>', '')
    cna_df['CNV_TYPE'] = cna_df.ALT_ALLELE.apply(__label_cnv_type)
    
    # Fill missing columns from mut_df
    cna_df.rename(columns={'CNV_CALL_ID' : 'VARIANT_CALL_ID',
                           'GENE' : 'CANONICAL_GENE'},
                  inplace=True)
    for prefix in 'CANONICAL HARMONIZED'.split():
        for suffix in 'PROTEIN_CHANGE CDNA_CHANGE VARIANT_CLASS'.split():
            cna_df['{}_{}'.format(prefix, suffix)] = cna_df.CNV_TYPE
    na_cols = 'CANONICAL_EXON CANONICAL_REF_SEQ_TSCP_ID CANONICAL_ENSEMBL_TSCP_ID ' + \
              'HARMONIZED_TRANSCRIPT_ID PATHOLOGIST_PATHOGENICITY ALLELE_FRACTION ' + \
              'COVERAGE MAX_GNOMAD_FREQUENCY'
    for colname in na_cols.split():
        cna_df[colname] = pd.NA
    cna_df.VARIANT_CALL_ID = 'CNV_' + cna_df.VARIANT_CALL_ID.astype(str)

    # Map approximate gene coordinates onto chrom/pos
    gvals = cna_df.CANONICAL_GENE.map(gene_map)
    backup_gvals = cna_df.HARMONIZED_HUGO_GENE_NAME.map(gene_map).to_dict()
    gvals.fillna(backup_gvals, inplace=True)
    gvals = gvals.apply(lambda v: v if isinstance(v, list) else [pd.NA, pd.NA])
    gvals_df = pd.DataFrame.from_dict(dict(zip(gvals.index, gvals.values)), 
                                      orient='index', dtype='object',
                                      columns='CHROMOSOME POSITION'.split())
    cna_df = cna_df.merge(gvals_df, how='left', left_index=True, right_index=True)

    # Final cleanup (subset to columns in mut_df, merge with mut_df, and resort by chrom/pos)
    out_df = pd.concat([cna_df.loc[:, mut_df.columns.tolist()], mut_df], ignore_index=True)

    return out_df.sort_values('CHROMOSOME POSITION'.split())


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutation-csv', help='Somatic mutation .csv', required=True)
    parser.add_argument('--cna-csv', help='Copy number alteration .csv', required=True)
    parser.add_argument('--samples-list', help='List of samples to evaluate', required=True)
    parser.add_argument('--id-map-tsv', help='.tsv mapping between various IDs', required=True)
    parser.add_argument('--genes-gtf', help='.gtf of gene annotations. If provided, ' +
                        'will be used to add approximate coordinates to CNA rows.')
    parser.add_argument('-o', '--outfile', default='stdout', help='path to somatic ' +
                        'variation .tsv [default: stdout]')
    args = parser.parse_args()

    # Load ID map for samples of interest
    id_map = load_id_map(args.samples_list, args.id_map_tsv)

    # Load mutation data and subset to samples of interest
    mut_df = load_mutation_data(args.mutation_csv, id_map)
    
    # Load CNA data, subset to samples of interest, and merge with mut_df
    mut_df = add_cna_data(args.cna_csv, mut_df, id_map, args.genes_gtf)

    # Write merged somatic data to output file
    if args.outfile in '- stdout'.split():
        fout = stdout
    else:
        fout = args.outfile
    mut_df = mut_df.rename(columns={'CHROMOSOME' : '#CHROMOSOME'})
    mut_df.to_csv(fout, sep='\t', index=False, na_rep='.')


if __name__ == '__main__':
    main()

