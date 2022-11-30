#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert COSMIC Cancer Mutation Census .tsv to VCF for VEP annotation
"""


import argparse
import pandas as pd
import pysam
import re
from somatic_df_to_vcf import format_alleles
from sys import stdout


chromsort_key = {str(i) : i  for i in range(1, 23, 1)}
chromsort_key['X'] = 23
chromsort_key['Y'] = 24


def load_header(header_in):
    """
    Load & extend input VCF header
    """

    # Load
    header = pysam.VariantFile(header_in).header
    # header.add_sample('dummy')

    # Add COSMIC INFO fields
    header.add_line('##INFO=<ID=COSMIC_GENE,Number=.,Type=String,Description="COSMIC gene">')
    header.add_line('##INFO=<ID=COSMIC_AA_CHANGE,Number=.,Type=String,Description="COSMIC amino acid change">')
    header.add_line('##INFO=<ID=COSMIC_GENE_TIER,Number=1,Type=Integer,Description="COSMIC gene tier">')
    header.add_line('##INFO=<ID=COSMIC_MUT_SIG,Number=1,Type=Integer,Description="COSMIC mutation significance tier">')
    header.add_line('##INFO=<ID=COSMIC_MUT_FREQ,Number=1,Type=Float,Description="COSMIC mutation frequency">')

    return header


def get_cmc_iterable(cmc_tsv):
    """
    Loads, preprocesses, and sorts CMC data
    Returns a generator of dicts, one per mutation
    """

    cols_to_keep = ['Mutation genome position GRCh37', 'GENE_NAME', 'CGC_TIER',
                    'COSMIC_SAMPLE_TESTED', 'COSMIC_SAMPLE_MUTATED', 
                    'MUTATION_SIGNIFICANCE_TIER', 'GENOMIC_WT_ALLELE_SEQ',
                    'GENOMIC_MUT_ALLELE_SEQ', 'Mutation AA']

    # Load CMC data
    df = pd.read_csv(cmc_tsv, sep='\t', usecols=cols_to_keep).\
            rename(columns={'GENOMIC_WT_ALLELE_SEQ' : 'REF_ALLELE',
                            'GENOMIC_MUT_ALLELE_SEQ' : 'ALT_ALLELE'})
    df.loc[:, 'REF_ALLELE ALT_ALLELE'.split()] \
        = df.loc[:, 'REF_ALLELE ALT_ALLELE'.split()].fillna('-')

    # Extract and format/correct coordinates
    coords = pd.DataFrame(df['Mutation genome position GRCh37'].\
                              astype(str).\
                              apply(lambda x: re.split(':|-', x)).\
                              tolist(), columns='CHROMOSOME POSITION END'.split())
    coords.CHROMOSOME = coords.CHROMOSOME.map({str(v) : k for k, v in chromsort_key.items()})
    coords['chr_idx'] = coords.CHROMOSOME.map(chromsort_key)
    out_df = pd.concat([coords, df], axis=1)
    out_df = out_df[~out_df.CHROMOSOME.isna() & ~out_df.POSITION.isna()]
    out_df.loc[:, 'POSITION END'.split()] \
        = out_df.loc[:, 'POSITION END'.split()].astype(int)
    ins_idx = (out_df.END - out_df.POSITION == 1)
    out_df.loc[ins_idx, 'POSITION'] = out_df.END[ins_idx]

    # Sort output values
    out_df = out_df.sort_values('chr_idx POSITION END'.split()).\
                    drop('chr_idx', axis=1)

    return out_df.iterrows()


def write_record(mvals, outvcf, ref_fa):
    """
    Format a mutation from CMC into a VCF record and write to outvcf
    """

    # Create record from basic variant info
    mvals = format_alleles(mvals, ref_fa, del_pos_adj=-1)
    new_id = '_'.join([mvals.CHROMOSOME, str(mvals.POSITION), *mvals.ALLELES])
    new_rec = outvcf.new_record(contig=mvals.CHROMOSOME, 
                                start=mvals.POSITION - 1, 
                                alleles=mvals.ALLELES, 
                                id=new_id)

    # Clean up INFO
    if len(new_rec.alleles[0]) == 1:
        if 'END' in new_rec.info.keys():
            new_rec.info.pop('END')
    if len(new_rec.alleles[1]) > len(new_rec.alleles[0]):
        new_rec.stop = mvals.END
    if not pd.isna(mvals.GENE_NAME):
        new_rec.info['COSMIC_GENE'] = mvals.GENE_NAME
    if not pd.isna(mvals['Mutation AA']):
        new_rec.info['COSMIC_AA_CHANGE'] = mvals['Mutation AA']
    if not pd.isna(mvals.CGC_TIER) \
    and mvals.CGC_TIER != 'Other':
        new_rec.info['COSMIC_GENE_TIER'] = int(mvals.CGC_TIER)
    if not pd.isna(mvals.MUTATION_SIGNIFICANCE_TIER) \
    and mvals.MUTATION_SIGNIFICANCE_TIER != 'Other':
        new_rec.info['COSMIC_MUT_SIG'] = int(mvals.MUTATION_SIGNIFICANCE_TIER)
    if not pd.isna(mvals.COSMIC_SAMPLE_MUTATED) \
    and not pd.isna(mvals.COSMIC_SAMPLE_TESTED) \
    and mvals.COSMIC_SAMPLE_TESTED > 0:
        new_rec.info['COSMIC_MUT_FREQ'] \
            = float(mvals.COSMIC_SAMPLE_MUTATED / mvals.COSMIC_SAMPLE_TESTED)

    # Write record to outvcf
    outvcf.write(new_rec)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cmc', help='COSMIC cmc_export.tsv.gz', required=True)
    parser.add_argument('--header', help='Header to use for output .vcf', required=True)
    parser.add_argument('--ref-fasta', help='Reference .fasta', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='output .vcf ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Load VCF header and open connection to output file
    header = load_header(args.header)

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.outfile, 'w', header=header)

    # Open connection to reference fa
    ref_fa = pysam.FastaFile(args.ref_fasta)

    # Iterate over sorted rows in CMC and write each to VCF
    for idx, mvals in get_cmc_iterable(args.cmc):
        write_record(mvals, outvcf, ref_fa)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

