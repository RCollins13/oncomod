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
import re


def mutdf_to_vcf(mut_df, samples, header_in, outfile, fasta_file, 
                 sample_column='DONOR_ID'):
    """
    Main wrapper function to convert mut_df to vcf
    """

    # Open connection to reference fa
    ref_fa = pysam.FastaFile(fasta_file)

    # Get list of samples to fill with no-call GTs (due to missing somatic variant data)
    no_mut_samples = samples.difference(set(mut_df[mut_df.END.isna()][sample_column].values))
    no_cna_samples = samples.difference(set(mut_df[~mut_df.END.isna()][sample_column].values))

    # Build dict of all records, metadata, and carrier samples
    vdict = compile_variant_dict(mut_df, ref_fa, sample_column)

    # Load header and add all sample IDs
    header = pysam.VariantFile(header_in).header
    for sample in samples:
        header.add_sample(sample)

    # Open connection to output VCF
    outvcf = pysam.VariantFile(outfile, 'w', header=header)

    # Iterate over unique variants and inject one VCF record for each into outvcf
    for vid, vdata in sorted(vdict.items(), key=lambda x: x[1]['ORDER']):
        # Skip variants on contigs not present in primary assembly
        if pd.isna(vdata['CHROM']) or pd.isna(vdata['POS']):
            continue
        if vdata['CHROM'] not in outvcf.header.contigs.keys():
            continue

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


def get_variant_id(series):
    """
    Build a unique variant ID from a single row in mut_df
    """

    vid_values = [re.sub('<|>', '', str(x)) for x in \
                  [series.CHROMOSOME, series.POSITION, *series.ALLELES]]
    
    return '_'.join(vid_values)


def format_alleles(item, ref_fa, del_pos_adj=0, verbose=False):
    """
    Get (ref, alt) allele tuple for a single variant
    Also correct ref/alt nomenclature for indels and update positions as needed
    """

    alleles = [item.REF_ALLELE, item.ALT_ALLELE.split(':')[0]]
    if alleles[1].startswith('<') and not alleles[1].endswith('>'):
        alleles[1] = alleles[1] + '>'
    alleles = tuple(alleles)

    # Handle insertions
    if alleles[0] == '-':
        new_ref = ref_fa.fetch(item.CHROMOSOME, item.POSITION, item.POSITION+1)
        new_alt = new_ref + alleles[1]
        alleles = (new_ref, new_alt)
        size = len(new_alt) - len(new_ref)
        if verbose:
            msg = 'Correcting insertion at {}:{}. Originally {}>{}, revised to {}>{}'
            print(msg.format(item.CHROMOSOME, item.POSITION, item.REF_ALLELE, 
                             item.ALT_ALLELE, new_ref, new_alt))

    # Handle deletions
    if alleles[1] == '-':
        new_ref = ref_fa.fetch(item.CHROMOSOME, item.POSITION-2, item.POSITION+len(alleles[0])-1)
        new_alt = new_ref[0]
        alleles = (new_ref, new_alt)
        size = len(new_ref) - len(new_alt)
        item['POSITION'] = item.POSITION + del_pos_adj
        item.END = item.POSITION + size
        if verbose:
            msg = 'Correcting deletion at {}:{}. Originally {}>{}, revised to {}>{}'
            print(msg.format(item.CHROMOSOME, item.POSITION, item.REF_ALLELE, 
                             item.ALT_ALLELE, new_ref, new_alt))

    item['ALLELES'] = alleles

    return item


def compile_variant_dict(mut_df, ref_fa, sample_column='DONOR_ID'):
    """
    Build a dict of all unique variants and their corresponding samples
    """

    variant_dict = {}
    i = 0

    for idx, item in mut_df.iterrows():
        
        item = format_alleles(item, ref_fa)
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
                                 'ALLELES' : item.ALLELES,
                                 'INFO' : {'SOURCE' : item.TEST_TYPE},
                                 'CARRIERS' : set([item[sample_column]])}
            if item.ALLELES[1] in '<DEL> <AMP>'.split():
                variant_dict[vid]['INFO']['SVTYPE'] = re.sub('<|>', '', item.ALLELES[1])
        else:
            variant_dict[vid]['CARRIERS'].add(item[sample_column])

    return variant_dict

