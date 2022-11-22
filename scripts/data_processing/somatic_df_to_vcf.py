#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Helper functions to convert a dataframe of somatic mutations to VCF

Note: this is not intended to be a standalone script. It is imported by 
preprocess_tcga_somatic.py and preprocess_dfci_profile_somatic.py. Refer to those
scripts for implementation details
"""


import pandas as pd
import pysam


def mutdf_to_vcf(mut_df, samples, header_in, outfile):
    """
    Main wrapper function to convert mut_df to vcf
    """

    # Get list of samples to fill with no-call GTs (due to missing somatic variant data)
    no_mut_samples = samples.difference(set(mut_df[mut_df.END.isna()].DONOR_ID.values))
    no_cna_samples = samples.difference(set(mut_df[~mut_df.END.isna()].DONOR_ID.values))

    # Build dict of all records, metadata, and carrier samples
    vdict = compile_variant_dict(mut_df)

    # Load header and add all sample IDs
    header = pysam.VariantFile(header_in).header
    for sample in samples:
        header.add_sample(sample)

    # Open connection to output VCF
    outvcf = pysam.VariantFile(outfile, 'w', header=header)

    # Iterate over unique variants and inject one VCF record for each into outvcf
    for vid, vdata in sorted(vdict.items(), key=lambda x: x[1]['ORDER']):
        new_rec = outvcf.new_record(contig=vdata['CHROM'], start=vdata['POS'], 
                                    stop=vdata['END'], alleles=vdata['ALLELES'], 
                                    id=vid)

        # Add INFO
        for key, value in vdata['INFO'].items():
            new_rec.info[key] = value

        # Add genotypes for each sample
        AN, AC = [0] * 2
        for sample in samples:
            # Leave samples missing mutation or CNA data as null
            if new_rec.alleles[1] in '<DEL> <AMP>'.split():
                if sample in no_cna_samples:
                    continue
            else:
                if sample in no_mut_samples:
                    continue

            # Assume samples are ref if they aren't carriers and aren't missing data
            AN += 1
            if sample in vdata['CARRIERS']:
                AC += 1
                new_rec.samples[sample]['GT'] = (1)
            else:
                new_rec.samples[sample]['GT'] = (0)

        # Add AC/AN/AF to INFO
        new_rec.info['AN'] = AN
        new_rec.info['AC'] = AC
        new_rec.info['AF'] = AC / AN

        # Write record to outvcf
        outvcf.write(new_rec)

    # Close outvcf to clear buffer
    outvcf.close()


def compile_variant_dict(mut_df):
    """
    Build a dict of all unique variants and their corresponding samples
    """

    variant_dict = {}
    i = 0

    for idx, item in mut_df.iterrows():
        
        vid = get_variant_id(item)

        if vid not in variant_dict.keys():
            i += 1
            if pd.isna(item.END):
                end_val = item.POSITION
            else:
                end_val = item.END
            variant_dict[vid] = {'ORDER' : i,
                                 'CHROM' : item.CHROMOSOME,
                                 'POS' : item.POSITION - 1,
                                 'END' : end_val,
                                 'ID' : vid,
                                 'ALLELES' : (item.REF_ALLELE, item.ALT_ALLELE),
                                 'INFO' : {'SOURCE' : item.TEST_TYPE},
                                 'CARRIERS' : set([item.DONOR_ID])}
        else:
            variant_dict[vid]['CARRIERS'].add(item.DONOR_ID)

    return variant_dict


def get_variant_id(series):
    """
    Build a unique variant ID from a single row in mut_df
    """

    if pd.isna(series.END):
        vid_cols = 'CHROMOSOME POSITION REF_ALLELE ALT_ALLELE'.split()
    else:
        vid_cols = 'CHROMOSOME POSITION END CANONICAL_CDNA_CHANGE'.split()

    return '_'.join(series[vid_cols].astype(str).values.tolist())

