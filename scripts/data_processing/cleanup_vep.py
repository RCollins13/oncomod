#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Clean up verbose VEP output for RASMod VCFs
"""


import argparse
import pandas as pd
import pysam
from copy import deepcopy
from sys import stdin, stdout


# Define values used in various functions below
vep_pop = 'Allele Existing_variation ALLELE_NUM STRAND UTRAnnotator_existing_InFrame_oORFs ' + \
          'UTRAnnotator_existing_OutOfFrame_oORFs UTRAnnotator_existing_uORFs ' + \
          'CADD_RAW LoF_filter LoF_flags LoF_info gnomADg gnomADe ClinVar COSMIC ' + \
          'IMPACT FLAGS MINIMISED SYMBOL_SOURCE HGVS_OFFSET SOURCE Feature'
vep_pop = vep_pop.split()
vep_severity = {'transcript_ablation' : 0,
                'splice_acceptor_variant' : 1,
                'splice_donor_variant' : 2,
                'stop_gained' : 3,
                'frameshift_variant' : 4,
                'stop_lost' : 5,
                'start_lost' : 6,
                'transcript_amplification' : 7,
                'inframe_insertion' : 8,
                'inframe_deletion' : 9,
                'missense_variant' : 10,
                'protein_altering_variant' : 11,
                'splice_region_variant' : 12,
                'splice_donor_5th_base_variant' : 13,
                'splice_donor_region_variant' : 14,
                'splice_polypyrimidine_tract_variant' : 15,
                'incomplete_terminal_codon_variant' : 16,
                'start_retained_variant' : 17,
                'stop_retained_variant' : 18,
                'synonymous_variant' : 19,
                'coding_sequence_variant' : 20,
                'mature_miRNA_variant' : 21,
                '5_prime_UTR_variant' : 22,
                '3_prime_UTR_variant' : 23,
                'non_coding_transcript_exon_variant' : 24,
                'intron_variant' : 25,
                'NMD_transcript_variant' : 26,
                'non_coding_transcript_variant' : 27,
                'upstream_gene_variant' : 28,
                'downstream_gene_variant' : 29,
                'TFBS_ablation' : 30,
                'TFBS_amplification' : 31,
                'TF_binding_site_variant' : 32,
                'regulatory_region_ablation' : 33,
                'regulatory_region_amplification' : 34,
                'feature_elongation' : 35,
                'regulatory_region_variant' : 36,
                'feature_truncation' : 37,
                'intergenic_variant' : 38}


def parse_vep_map(invcf):
    """
    Parse VEP field mappings to variable names
    """

    vep_text = invcf.header.info.get('CSQ').description
    vep_fields = vep_text.split('Format: ')[1].replace('\'', '').split('|')

    return {i : k for i, k in enumerate(vep_fields)}


def reformat_header(invcf):
    """
    Reformat header for output VCF
    """

    out_header = invcf.header

    # Modify VEP CSQ entry to reflect removed values
    old_vep = out_header.info.get('CSQ')
    old_descr = old_vep.description
    old_fields = old_descr.split('Format: ')[1].split('|')
    new_fields = [f for f in old_fields if f not in vep_pop]
    new_descr = old_descr.split('Format: ')[0] + 'Format: ' + '|'.join(new_fields)
    out_header.info.remove_header('CSQ')
    out_header.add_meta(key='INFO', items=[('ID', 'CSQ'), ('Number', '.'),
                                           ('Type', 'String'), 
                                           ('Description', new_descr)])

    return out_header


def load_tx_map(tx_tsv):
    """
    Load --transcript-info as a dict
    """

    import pdb; pdb.set_trace()
    
    tx_df = pd.DataFrame(tx_tsv, sep='\t')
    
    return {'ENSG' : ,
            'symbol' : ,
            'tx_len' : }


def cleanup(record, vep_map, tx_map):
    """
    Simplify VEP entry for a single record
    """

    # Build pd.DataFrame of all VEP entries
    # (Excluding keys in vep_pop, defined above)
    vep_vals = {}
    for i, vep_str in enumerate(record.info.get('CSQ')):
        vep_vals[i] = {k : v for k, v in zip(vep_map.values(), vep_str.split('|')) \
                       if k not in vep_pop}
    vdf = pd.DataFrame.from_dict(vep_vals, orient='index')

    # Select single VEP entry to retain per gene per variant
    # Priority:
    #   1. Always choose protein-coding over non-coding transcripts
    #   2. Take most severe consequence, if multiple are present
    #   3a. Take canonical transcript, if most severe consequence present on multiple
    #   3b. Take longest transcript, if most severe consequence present on multiple
    keep_idxs = set()
    for gene in set(vdf.Gene.values):

        # Check for protein-coding transcripts
        coding = (vdf.BIOTYPE == 'protein_coding')
        # if any(coding):

        # vdf[vdf.Gene == gene]
        import pdb; pdb.set_trace()



def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input VEP-annotated .vcf [default: stdin]',
                        default='stdin')
    parser.add_argument('vcf_out', help='output .vcf [default: stdout]', default='stdout')
    parser.add_argument('-t', '--transcript-info', required=True, help='.tsv ' + \
                        'mapping ENST:ENSG:symbol:length')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Parse mapping of VEP fields and transcript info
    vep_map = parse_vep_map(invcf)
    tx_map = load_tx_map(args.transcript_info)
    
    # Open connection to output file
    out_header = reformat_header(invcf)
    if args.vcf_out in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=out_header)

    # Iterate over records in invcf, clean up each, and write to outvcf
    for record in invcf.fetch():
        record = cleanup(record, vep_map, tx_map)
        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

