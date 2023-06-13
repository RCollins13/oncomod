#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Merge OncoPanel mutation and CNA data to generate input for All-FIT
"""


import argparse
import pandas as pd


autosomes = []
for i in range(22):
    autosomes.append(str(i+1))
    autosomes.append('chr{}'.format(i+1))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutations', help='OncDRS GENOMIC_MUTATION_RESULTS.csv',
                        required=True)
    parser.add_argument('--cnas', help='OncDRS GENOMIC_CNV_RESULTS.csv', 
                        required=True)
    parser.add_argument('-o', '--outdir', help='path to output directory',
                        default='./')
    args = parser.parse_args()

    # Load mutations and retain relevant columns
    # Note that All-FIT only works on autosomes, so a chromosome filter is applied
    mut_cols = 'SAMPLE_ACCESSION_NBR VARIANT_CALL_ID CANONICAL_GENE' + \
               ' COVERAGE ALLELE_FRACTION CHROMOSOME'
    mut = pd.read_csv(args.mutations, sep=',', usecols=mut_cols.split())
    mut = mut[mut.CHROMOSOME.astype(str).str.rstrip().isin(autosomes)]
    mut.drop(columns='CHROMOSOME', inplace=True)
    mut.rename(columns = {'CANONICAL_GENE' : 'GENE'}, inplace=True)

    # Load CNAs and retain the relevant columns
    cna_cols = 'SAMPLE_ACCESSION_NBR GENE COPY_COUNT'.split()
    cna = pd.read_csv(args.cnas, sep=',', usecols=cna_cols)
    
    # Annotate ploidy for all mutations based on matching sample accession & gene
    # Fill missing copy number information with naive diploid prior
    # Round negative estimated copy numbers to zero
    out_df = mut.merge(cna, on='SAMPLE_ACCESSION_NBR GENE'.split(), how='left')
    out_df.COPY_COUNT = out_df.COPY_COUNT.fillna(value=2).astype('int')
    out_df.loc[out_df.COPY_COUNT < 0, 'COPY_COUNT'] = 0

    # Clean data by excluding samples missing AF or depth info
    out_df = out_df.loc[(~out_df.COPY_COUNT.isna()) & (~out_df.COVERAGE.isna()), :]
    out_df.rename(columns={'VARIANT_CALL_ID' : 'ID', 'ALLELE_FRACTION' : 'Allele_Freq',
                           'COVERAGE' : 'Depth', 'COPY_COUNT' : 'Ploidy'},
                  inplace=True)
    
    # Write one All-FIT input file per accession number
    out_cols = 'ID Allele_Freq Depth Ploidy'.split()
    sids = out_df.SAMPLE_ACCESSION_NBR.unique().tolist()
    print('Identified at least one variant for {:,} samples'.format(len(sids)))
    log_fmt = 'Processed {:,} of {:,} samples ({:.1f}%)'
    k = 0
    for sid in sids:
        k += 1
        out_path = args.outdir + '/' + sid + '.AllFIT_input.tsv'
        keepers = out_df.SAMPLE_ACCESSION_NBR == sid
        if len(keepers) > 0:
            out_df.loc[keepers, out_cols].to_csv(out_path, index=False, sep='\t')
        if k % 100 == 0:
            print(log_fmt.format(k, len(sids), 100 * k / len(sids)))


if __name__ == '__main__':
    main()

