#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert GTEx eQTL gene-variant pair data to VCF for VEP annotation
"""


import argparse
import pandas as pd
import pysam
import re
from os.path import basename
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

    # Add GTEx INFO fields
    header.add_line('##INFO=<ID=GTEx_eGene,Number=.,Type=String,Description="List of ENSG IDs for GTEx genes whose expression is significantly predicted by this variant.">')
    header.add_line('##INFO=<ID=GTEx_eQTL_beta,Number=.,Type=Float,Description="List of expression effect size estimates in GTEx. One per entry in GTEx_eGene.">')
    header.add_line('##INFO=<ID=GTEx_eQTL_tissue,Number=.,Type=String,Description="List of tissues with significant eQTL effects in GTEx. One per entry in GTEx_eGene.">')

    return header


def add_tissue(var_dict, tsv_in):
    """
    Add eQTL results from a single tissue to dictionary of all variants
    """

    # Read data
    keep_cols = 'variant_id gene_id slope'
    new_df = pd.read_csv(tsv_in, sep='\t', usecols=keep_cols.split())

    # Infer tissue from filename
    new_df['tissue'] = basename(tsv_in).split('.')[0]

    # Process each row in new_df
    for idx, vdat in new_df.iterrows():

        # Get basic variant info and add to var_dict if not already present
        vid_parts = vdat.variant_id.split('_')
        vid = '_'.join(vid_parts[:4])
        if vid not in var_dict.keys():
            chrom, pos, ref, alt = vid_parts[:4]
            var_dict[vid] = {'chrom' : chrom, 'pos' : pos, 'ref' : ref, 'alt' : alt,
                             'eGenes' : [], 'betas' : [], 'tissues' : []}

        # Add eQTL info to var_dict
        var_dict[vid]['eGenes'].append(vdat.gene_id)
        var_dict[vid]['betas'].append(vdat.slope)
        var_dict[vid]['tissues'].append(vdat.tissue)

    return var_dict


def write_record(vdat, vid, outvcf):
    """
    Format a GTEx eQTL variant into a VCF record and write to outvcf
    """

    # Create record from basic variant info
    new_rec = outvcf.new_record(contig=vdat['chrom'], 
                                start=int(vdat['pos']) - 1, 
                                alleles=(vdat['ref'], vdat['alt']), 
                                id=vid)

    # Write eQTL data to INFO
    new_rec.info['GTEx_eGene'] = vdat['eGenes']
    new_rec.info['GTEx_eQTL_beta'] = vdat['betas']
    new_rec.info['GTEx_eQTL_tissue'] = vdat['tissues']

    # Write record to outvcf
    outvcf.write(new_rec)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtex-tsv', required=True, action='append',
                        help='GTEx gene-variant .tsv. May be specified multiple times.')
    parser.add_argument('--header', help='Header to use for output .vcf', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='output .vcf ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Load VCF header
    header = load_header(args.header)

    # Build dict of all variants with collapsed information for each
    var_dict = {}
    for tsv_in in args.gtex_tsv:
        var_dict = add_tissue(var_dict, tsv_in)

    # Open connection to output file
    if args.outfile in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.outfile, 'w', header=header)

    # Iterate over sorted entries in var_dict and write each to VCF
    for vid, vdat in var_dict.items():
        write_record(vdat, vid, outvcf)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

