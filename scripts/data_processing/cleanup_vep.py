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
import numpy as np
import pandas as pd
import pysam
from copy import deepcopy
from sys import stdin, stdout


# Define values used in various functions below
vep_pop = 'Allele Existing_variation ALLELE_NUM STRAND UTRAnnotator_existing_InFrame_oORFs ' + \
          'UTRAnnotator_existing_OutOfFrame_oORFs UTRAnnotator_existing_uORFs ' + \
          'CADD_RAW LoF_filter LoF_flags LoF_info gnomADg gnomADe ClinVar COSMIC ' + \
          'IMPACT FLAGS MINIMISED SYMBOL_SOURCE HGVS_OFFSET SOURCE'
vep_pop = vep_pop.split()
vep_reloc_map = {'CADD_PHRED' : 'CADD',
                 'GERP' : 'GERP',
                 'phastCons' : 'phastCons',
                 'phyloP' : 'phyloP',
                 'COSMIC_COSMIC_MUT_FREQ' : 'COSMIC_freq'}
vep_reloc = list(vep_reloc_map.keys())
gnomad_pops = 'AFR AMR ASJ EAS FIN NFE OTH POPMAX'.split()
ensembl_reg = 'MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE ' + \
              'TRANSCRIPTION_FACTORS CELL_TYPE'
ensembl_reg = ensembl_reg.split()
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
clinvar_sig = {'Pathogenic' : 0,
               'Pathogenic/Likely_pathogenic' : 1,
               'Likely_pathogenic' : 2,
               'Uncertain_significance' : 3,
               'Conflicting_interpretations_of_pathogenicity' : 4,
               'Likely_benign' : 5,
               'Benign' : 6}
clinvar_remap = {'Pathogenic' : 'P',
                'Pathogenic/Likely_pathogenic' : 'PLP',
                'Likely_pathogenic' : 'LP',
                'Uncertain_significance' : 'UNCERTAIN',
                'Conflicting_interpretations_of_pathogenicity' : 'CONFLICTING',
                'Likely_benign' : 'LB',
                'Benign' : 'B'}


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
    cols_to_drop = vep_pop + vep_reloc + ensembl_reg
    for pop in gnomad_pops:
        cols_to_drop += 'gnomADe_AF_{} gnomADg_AF_{}'.format(pop, pop).split()
    old_vep = out_header.info.get('CSQ')
    old_descr = old_vep.description
    old_fields = old_descr.split('Format: ')[1].split('|')
    new_fields = [f for f in old_fields if f not in cols_to_drop]
    new_descr = old_descr.split('Format: ')[0] + 'Format: ' + '|'.join(new_fields)
    out_header.info.remove_header('CSQ')
    out_header.add_meta(key='INFO', items=[('ID', 'CSQ'), ('Number', '.'),
                                           ('Type', 'String'), 
                                           ('Description', new_descr)])

    # Add new header lines for values relocated from CSQ
    for key in vep_reloc_map.values():
        if key in out_header.info.keys():
            out_header.info.remove_header(key)
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'CADD'), ('Number', '1'), ('Type', 'Float'),
                               ('Description', 'Phred-scaled CADD score')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'GERP'), ('Number', '1'), ('Type', 'Float'),
                               ('Description', 'GERP score')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'phastCons'), ('Number', '1'), ('Type', 'Float'),
                               ('Description', 'phastCons score')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'phyloP'), ('Number', '1'), ('Type', 'Float'),
                               ('Description', 'phyloP score')])
    out_header.info.remove_header('ClinVar')
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'ClinVar'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'Highest ClinVar significance ' + \
                                'designation')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'COSMIC_mutation_tier'), ('Number', '1'), ('Type', 'Integer'),
                               ('Description', 'Highest COSMIC mutation tier')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'COSMIC_freq'), ('Number', '1'), ('Type', 'Float'),
                               ('Description', 'Pan-cancer mutation frequency reported ' + \
                                'by COSMIC')])

    # Clean up gnomAD-related fields
    for suffix in 'e g'.split():
        out_header.info.remove_header('gnomAD' + suffix)
        for pop in gnomad_pops:
            out_header.info.remove_header('gnomAD{}_AF_{}'.format(suffix, pop))
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'gnomAD_source'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'gnomAD allele frequency source (exomes or genomes)')])
    for pop in gnomad_pops:
        if pop != 'POPMAX':
            descr = 'gnomAD allele frequency in {} samples in gnomAD v2.0.1'.format(pop)
        else:
            descr = 'Maximum allele frequency among all populations in gnomAD v2.0.1'
        out_header.add_meta(key='INFO', 
                            items=[('ID', 'gnomAD_AF_' + pop), ('Number', '1'), ('Type', 'Float'),
                                   ('Description', descr)])

    return out_header


def load_tx_map(tx_tsv):
    """
    Load --transcript-info as a dict
    """
    
    # Load data
    tx_df = pd.read_csv(tx_tsv, sep='\t', header=None)
    tx_df.columns = 'ENST ENSG symbol tx_len'.split()

    # Remove Ensembl version info from ENSG and ENST IDs
    for ecol in 'ENST ENSG'.split():
        tx_df[ecol] = tx_df[ecol].str.replace('.[0-9]+$', '', regex=True)

    # Map everything vs. ENST IDs
    tx_df.set_index('ENST', drop=True, inplace=True)
    
    return {'ENSG' : tx_df.ENSG.to_dict(),
            'ENSG_to_symbol' : tx_df.symbol.set_axis(tx_df.ENSG).to_dict(),
            'symbol' : tx_df.symbol.to_dict(),
            'tx_len' : tx_df.tx_len.to_dict()}


def _harmonize_gnomad(record, vdf):
    """
    Harmonize allele frequencies between gnomAD exomes and genomes
    Use exome data where possible due to larger sample size
    """

    exome_cols = [k for k in vdf.columns if 'gnomADe_' in k]
    genome_cols = [k for k in vdf.columns if 'gnomADg_' in k]
    has_exome = any((vdf[exome_cols] != '').transpose())
    has_genome = any((vdf[genome_cols] != '').transpose())
    if has_exome:
        gsrc = 'exome'
    elif has_genome:
        gsrc = 'genome'
    else:
        gsrc = 'neither'
    record.info['gnomAD_source'] = gsrc
    for pop in gnomad_pops:
        if gsrc != 'neither':
            af_vals = vdf['gnomAD{}_AF_{}'.format(gsrc[0], pop)]
            af = af_vals[~af_vals.isin(['', '.'])].astype(float).mean()
            if pd.isna(af):
                af = 0
        else:
            af = 0
        try:
            record.info['gnomAD_AF_' + pop] = af
        except:
            import pdb; pdb.set_trace()

    vdf.drop(exome_cols + genome_cols, axis=1, inplace=True)

    return record, vdf


def _clean_regulatory(record, vdf, tx_map):
    """
    Move Ensembl regulatory annotations from VEP CSQ field into INFO
    """

    # Process each regulatory entry separately based on BIOTYPE
    idx_to_drop = set()
    for idx, row in vdf[vdf.Consequence == 'regulatory_region_variant'].iterrows():

        idx_to_drop.add(idx)

        # Skip promoters (these are annotated using a custom BED file)
        promoter_biotypes = 'promoter promoter_flanking_region'
        if row.BIOTYPE in promoter_biotypes.split():
            continue

        # Check if feature is annotated in any tissues of interest
        # ts_biotypes = 'TF_binding_site CTCF_binding_site open_chromatin_region enhancer'
        # if row.BIOTYPE in ts_biotypes.split():
        activity = {x.split(':')[0] : x.split(':')[1] for x in row.CELL_TYPE.split('&')}
        if not all(pd.Series(activity.values()).isin('NA INACTIVE'.split())):
            import pdb; pdb.set_trace()

        # else:
        #     import pdb; pdb.set_trace()

    # Remove Ensembl regulatory entries from CSQ
    vdf.drop(ensembl_reg, axis=1, inplace=True)
    vdf.drop(idx_to_drop, axis=0, inplace=True)

    return record, vdf


def cleanup(record, vep_map, tx_map):
    """
    Simplify VEP entry for a single record
    """

    # Build pd.DataFrame of all VEP entries
    # (Excluding keys in vep_pop, defined above)
    vep_vals = {}
    for i, vep_str in enumerate(record.info.get('CSQ', [])):
        vep_vals[i] = {k : v for k, v in zip(vep_map.values(), vep_str.split('|')) \
                       if k not in vep_pop}
    vdf = pd.DataFrame.from_dict(vep_vals, orient='index')

    # Don't process records lacking CSQ INFO field
    if len(vdf) == 0:
        return record

    # Pre-filter VEP entries to ignore all intergenic entries if any genic entries are present
    genic_csqs = [k for k, v in vep_severity.items() if v <= 27]
    genic_hits = vdf.Consequence.isin(genic_csqs)
    intergenic_csqs = [k for k, v in vep_severity.items() if v > 27]
    intergenic_hits = vdf.Consequence.isin(intergenic_csqs)
    if any(genic_hits) and any(intergenic_hits):
        vdf = vdf[genic_hits]

    # Select single VEP entry to retain per gene per variant
    # Priority:
    #   1. Always choose protein-coding over non-coding transcripts
    #   2. Take most severe consequence, if multiple are present
    #   3a. Take canonical transcript, if most severe consequence present on multiple
    #   3b. Take longest transcript, if most severe consequence present on multiple
    keep_idxs = set()
    for gene in set(vdf.Gene.values):

        gdf = vdf.loc[vdf.Gene == gene]

        # Subset to protein-coding transcripts, if multiple entries are present
        coding = (gdf.BIOTYPE == 'protein_coding')
        if any(coding) and len(gdf) > 1:
            gdf = gdf.loc[coding]

        # Subset to most severe consequence, if multiple are present
        gdf.loc[:, 'vep_priority'] = gdf.loc[:, 'Consequence'].map(vep_severity)
        if len(gdf.vep_priority.unique()) > 1:
            gdf = gdf[gdf.vep_priority == np.nanmin(gdf['vep_priority'])]

        # Check if any canonical transcripts are present
        canon = (gdf.CANONICAL == 'YES')
        if any(canon) and len(gdf) > 1:
            gdf = gdf.loc[canon]

        # Use transcript length as the final tiebreaker
        if len(gdf) > 1:
            gdf['tx_len'] = gdf.Feature.map(tx_map['tx_len'])
            gdf = gdf.sort_values('tx_len', ascending=False)

        # Write index to keep to keep_idxs
        keep_idxs.add(gdf.head(1).index[0])
    
    # Subset VEP entries to only those retained by the above criteria
    vdf = vdf.loc[keep_idxs, :]

    # Relocate values defined at a site level to INFO (and out of CSQ)
    for oldkey, newkey in vep_reloc_map.items():
        oldvals = vdf.loc[vdf[oldkey] != '', oldkey]
        if any(oldvals.str.contains('&')):
            oldvals = oldvals.str.split('&', expand=True).max(axis=1)
        newval = oldvals.astype(float).mean()
        if not pd.isna(newval):
            record.info[newkey] = newval
    vdf.drop(vep_reloc_map.keys(), axis=1, inplace=True)

    # Harmonize gnomAD frequencies between exomes and genomes
    record, vdf = _harmonize_gnomad(record, vdf)

    # Extract top ClinVar significance, if any reported
    if any(vdf.ClinVar_CLNSIG != ''):
        best_clinvar_idx = vdf.ClinVar_CLNSIG.map(clinvar_sig).sort_values().index.values[0]
        best_clinvar_label = clinvar_remap.get(vdf.ClinVar_CLNSIG[best_clinvar_idx])
        record.info['ClinVar'] = best_clinvar_label

    # Extract top COSMIC significance, if any reported
    if any(vdf.COSMIC_COSMIC_MUT_SIG != ''):
        best_cosmic = int(vdf.COSMIC_COSMIC_MUT_SIG.sort_values().head(1).values[0])
        record.info['COSMIC_mutation_tier'] = best_cosmic

    # Clean up Ensembl regulatory annotations
    record, vdf = _clean_regulatory(record, vdf, tx_map)

    # Overwrite CSQ INFO field with cleaned VEP info
    record.info.pop('CSQ')
    record.info['CSQ'] = tuple(['|'.join(vals) for i, vals in vdf.iterrows()])
    
    # Return cleaned record
    return record


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

