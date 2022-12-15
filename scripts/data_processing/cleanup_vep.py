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
import os
import pandas as pd
import pybedtools as pbt
import pysam
from copy import deepcopy
from sys import stdin, stdout, stderr, path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from general_utils import load_tx_map
from vep_utils import parse_vep_map, vep2df


# Define values used in various functions below
vep_pop = 'Allele Existing_variation ALLELE_NUM STRAND MOTIF_NAME MOTIF_POS ' + \
          'HIGH_INF_POS MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS ' + \
          'UTRAnnotator_existing_InFrame_oORFs ' + \
          'UTRAnnotator_existing_OutOfFrame_oORFs UTRAnnotator_existing_uORFs ' + \
          'CADD_RAW LoF_filter LoF_flags LoF_info gnomADg gnomADe ClinVar COSMIC ' + \
          'IMPACT FLAGS MINIMISED SYMBOL_SOURCE HGVS_OFFSET SOURCE'
vep_pop = vep_pop.split()
vep_reloc_map = {'CADD_PHRED' : 'CADD',
                 'GERP' : 'GERP',
                 'phastCons' : 'phastCons',
                 'phyloP' : 'phyloP',
                 'COSMIC_COSMIC_MUT_FREQ' : 'COSMIC_freq'}
vep_reloc = list(vep_reloc_map.keys()) + 'promoters ABC_enhancers'.split()
gnomad_pops = 'AFR AMR ASJ EAS FIN NFE OTH POPMAX'.split()
ensembl_reg = 'MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE ' + \
              'TRANSCRIPTION_FACTORS CELL_TYPE'
ensembl_reg = ensembl_reg.split()
ensembl_tissue_map = {'pancreas' : 'pancreas',
                      'endocrine_pancreas' : 'pancreas',
                      'sigmoid_colon' : 'colon',
                      'foreskin_melanocyte_1' : 'skin',
                      'foreskin_melanocyte_2' : 'skin',
                      'lung_1' : 'lung',
                      'lung_2' : 'lung'}
ensembl_reg_activity = 'ACTIVE POISED REPRESSED'.split()
spliceai_fields = 'SpliceAI_pred_DP_AG SpliceAI_pred_DP_AL SpliceAI_pred_DP_DG ' + \
                  'SpliceAI_pred_DP_DL SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL ' + \
                  'SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL SpliceAI_pred_SYMBOL'
spliceai_fields = spliceai_fields.split()
gtex_fields = 'GTEx GTEx_GTEx_eGene GTEx_GTEx_eQTL_beta GTEx_GTEx_eQTL_tissue'.split()
misc_header_info_to_pop = 'CCR ExAC_regional_constraint'.split()
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


def reformat_header(invcf):
    """
    Reformat header for output VCF
    """

    out_header = invcf.header

    # Modify VEP CSQ entry to reflect removed values
    cols_to_drop = vep_pop + vep_reloc + ensembl_reg + spliceai_fields + gtex_fields
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
    for key in vep_reloc + gtex_fields + misc_header_info_to_pop:
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
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'intergenic'), ('Number', '0'), ('Type', 'Flag'),
                               ('Description', 'Variant is intergenic.')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'RAS_distance'), ('Number', '1'), ('Type', 'Integer'),
                               ('Description', 'Distance to nearest RAS gene (for ' + \
                                'intergenic variants only).')])
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
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'promoter'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'List of genes whose promoter is ' + \
                                'overlapped by this variant')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'enhancer'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'List of all genes predicted to ' + \
                                'be regulated by an enhancer overlapping this ' + \
                                'variant. Format: gene|ABC_score|tissue')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'SpliceAI'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'List of all genes with a splicing ' + \
                                'consequence predicted by SpliceAI. Format: ' + \
                                'gene:pipe-delimited consequences')])
    out_header.add_meta(key='INFO', 
                        items=[('ID', 'GTEx_eQTL'), ('Number', '.'), ('Type', 'String'),
                               ('Description', 'List of all eQTL associations ' + \
                                'reported by GTEx for this variant. Format: ' + \
                                'gene|tissue|beta')])
    if 'CELL_TYPE' in old_fields:
        out_header.add_meta(key='INFO', 
                            items=[('ID', 'Ensembl_regulatory_feature'), 
                                   ('Number', '.'), ('Type', 'String'),
                                   ('Description', 'Comma-delimited list of feature type(s) ' + \
                                    'present at this locus in Ensembl regulatory build v1.0')])
        activity_descr_template = 'Comma-delimited list of cell types with {} marks ' + \
                                  'in Ensembl regulatory build v1.0'
        for activity in ensembl_reg_activity:
            out_header.add_meta(key='INFO', 
                                items=[('ID', 'Ensembl_{}_celltypes'.format(activity.lower())), 
                                       ('Number', '.'), ('Type', 'String'),
                                       ('Description', activity_descr_template.format(activity.lower()))])

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


def load_ras_bt(gtf_in):
    """
    Extract start/stop coordinates for KRAS/NRAS/HRAS
    """

    return pbt.BedTool(gtf_in).\
               filter(lambda x: x[2] == 'gene').\
               filter(lambda x: x.attrs.get('gene_name') in 'HRAS NRAS KRAS'.split()).\
               saveas()


def _harmonize_gnomad(record, vdf):
    """
    Harmonize allele frequencies between gnomAD exomes and genomes
    Use exome data where possible due to larger sample size
    """

    exome_cols = [k for k in vdf.columns if 'gnomADe_' in k]
    genome_cols = [k for k in vdf.columns if 'gnomADg_' in k]
    has_exome = (vdf[exome_cols] != '').any(axis=1).values[0]
    has_genome = (vdf[genome_cols] != '').any(axis=1).values[0]
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
        record.info['gnomAD_AF_' + pop] = af

    vdf.drop(exome_cols + genome_cols, axis=1, inplace=True)

    return record, vdf


def _clean_abc(record, vdf):
    """
    Relocate ABC enhancer data from INFO:CSQ to INFO
    """

    # Extract ABC enhancer-gene links
    abc_str = '&'.join(vdf.ABC_enhancers[~vdf.ABC_enhancers.duplicated()].values)
    abc_df = pd.Series(list(set(abc_str.split('&')))).str.split(':', expand=True)
    abc_df.columns = 'gene score sample tissue'.split()

    # For genes with enhancers in multiple samples from the same tissue, report 
    # the maximum ABC score
    abc_df = abc_df.drop('sample', axis=1).sort_values('score', ascending=False)
    abc_df = abc_df[~abc_df['gene tissue'.split()].duplicated()]

    # Reformat for single entry in INFO
    newdat = ','.join(['|'.join(v.values) for i, v in abc_df.iterrows()])

    # Update record and drop columns from vdf
    record.info['enhancer'] = newdat
    vdf.drop('ABC_enhancers', axis=1, inplace=True)

    return record, vdf


def _clean_spliceai(record, vdf, cutoff=0.5):
    """
    Relocalize SpliceAI predictions from INFO:CSQ to INFO
    """

    spliceai_pred_cols = vdf.columns[vdf.columns.str.contains('SpliceAI_pred_DS_')]
    all_splice_scores = vdf.loc[:, spliceai_pred_cols]
    all_splice_scores = all_splice_scores[~all_splice_scores.apply(lambda x: any(x == ''), axis=1)]
    splice_hits = all_splice_scores.astype(float) >= cutoff
    all_splice_csqs = all_splice_scores.columns.str.replace('SpliceAI_pred_DS_', '')
    splice_csqs = ['|'.join(all_splice_csqs[hits]) for i, hits in splice_hits.iterrows()]
    splice_info = ','.join([':'.join(v) for v in 
                            zip(vdf.SpliceAI_pred_SYMBOL, splice_csqs) 
                            if v[1] != ''])
    if splice_info != '':
        record.info['SpliceAI'] = splice_info
    vdf.drop([c for c in vdf.columns if 'SpliceAI' in c], inplace=True, axis=1)

    return record, vdf


def _clean_gtex(record, vdf, tx_map):
    """
    Relocate GTEx eQTL data from INFO:CSQ to INFO
    """

    gtex_tissue_map = {'Pancreas' : 'pancreas',
                       'Skin_Sun_Exposed_Lower_leg' : 'skin',
                       'Skin_Not_Sun_Exposed_Suprapubic' : 'skin',
                       'Colon_Sigmoid' : 'colon',
                       'Colon_Transverse' : 'colon',
                       'Lung' : 'lung'}

    # Extract GTEx eQTL data
    gtex_dat = vdf[gtex_fields].drop('GTEx', axis=1).head(1).values[0]
    ensgs, betas, tissues = [v.split('&') for v in gtex_dat]
    tissues = [gtex_tissue_map[t] for t in tissues]
    genes = [tx_map['ENSG_to_symbol'].get(e.split('.')[0]) for e in ensgs]
    gtex_df = pd.DataFrame([genes, tissues, betas], index='gene tissue beta'.split()).transpose()

    # For genes with effects in multiple samples from the same tissue, report the 
    # absolute maximum effect size
    gtex_df['abs_beta'] = gtex_df.beta.astype(float).abs()
    gtex_df.sort_values('abs_beta', ascending=False, inplace=True)
    gtex_df = gtex_df[~gtex_df['gene tissue'.split()].duplicated()].drop('abs_beta', axis=1)

    # Reformat for single entry in INFO
    newdat = ','.join(['|'.join(v.values) for i, v in gtex_df.iterrows()])

    # Update record and drop columns from vdf
    record.info['GTEx_eQTL'] = newdat
    vdf.drop(gtex_fields, axis=1, inplace=True)

    return record, vdf


def _update_str_set(record, key, newvals, sep=','):
    """
    Updates an INFO field in a VCF record by adding one or more newvals to the existing values
    Will add directly if the field is not already present
    """

    if key not in record.info.keys():
        record.info[key] = sep.join(set(newvals))
    else:
        oldvals = set(record.info[key].split(sep))
        oldvals.add(newvals)
        record.info[key] = sep.join(oldvals)

    return record


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
        activity = {x.split(':')[0] : x.split(':')[1] for x in row.CELL_TYPE.split('&')}
        if any(pd.Series(activity.values()).isin(ensembl_reg_activity)):

            # Add regulatory feature to INFO field
            record = _update_str_set(record, 'Ensembl_regulatory_feature', [row.BIOTYPE])

            # Annotate active/poised/repressed cell types
            act_fmt = 'Ensembl_{}_celltypes'
            for ct, act in activity.items():
                if act not in ensembl_reg_activity:
                    continue
                record = _update_str_set(record, act_fmt.format(act.lower()), [ct])

    # Remove Ensembl regulatory entries from CSQ
    vdf.drop(set(ensembl_reg).intersection(set(vdf.columns)), 
             axis=1, inplace=True)
    vdf.drop(idx_to_drop, axis=0, inplace=True)

    return record, vdf


def cleanup_cnv(record, vep_fields, gtf_tabix):
    """
    Add simple gene overlap CSQ entry to CNVs
    """

    if 'chr' in record.chrom:
        chrom = record.chrom
    else:
        chrom = 'chr' + record.chrom

    # Get all annotated genes overlapping CNV
    hits = [g.rstrip().split('\t') for g in \
            gtf_tabix.fetch(chrom, record.pos, record.stop) \
            if g.split('\t')[2] == 'gene']
    
    # Extract GTF annotations for each as dict
    hit_dicts = []
    for hit in hits:
        hit_dicts.append({v.split()[0] : v.split()[1] for v in 
                          hit[8].replace('"', '').split('; ')})

    # Build CSQ entries for each hit
    csq_entries = []
    for hit_dict in hit_dicts:
        vep_vals = {k : '' for k in vep_fields}

        # Set CNV type-agnostic values
        vep_vals['SYMBOL'] = hit_dict.get('gene_name', '')
        vep_vals['Gene'] = hit_dict.get('gene_id', '').split('.')[0]
        vep_vals['Feature_type'] = 'Transcript'
        vep_vals['Feature'] = hit_dict.get('transcript_id', '').split('.')[0]
        vep_vals['BIOTYPE'] = hit_dict.get('transcript_type', '')

        # Set CNV type-specific values
        if record.info['SVTYPE'] in 'AMP DUP'.split():
            vep_vals['Consequence'] = 'transcript_amplification'

        elif record.info['SVTYPE'] == 'DEL':
            vep_vals['Consequence'] = 'transcript_ablation'

        csq_entries.append('|'.join(vep_vals.values()))

    record.info['CSQ'] = ','.join(csq_entries)

    return record


def cleanup(record, vep_map, tx_map, ras_bt, spliceai_cutoff=0.5):
    """
    Simplify VEP entry for a single record
    """

    # Build pd.DataFrame of all VEP entries
    # (Excluding keys in vep_pop, defined above)
    vdf = vep2df(record, vep_map, vep_pop)

    # Don't process records lacking CSQ INFO field
    if len(vdf) == 0:
        return record

    # Pre-filter VEP entries to ignore all intergenic entries if any genic entries are present
    genic_csqs = [k for k, v in vep_severity.items() if v <= 27]
    genic_hits = vdf.Consequence.isin(genic_csqs)
    intergenic_csqs = [k for k, v in vep_severity.items() if v > 27]
    intergenic_hits = vdf.Consequence.isin(intergenic_csqs)
    if genic_hits.any():
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

    # Annotate distance to nearest RAS gene for purely intergenic variants
    if vdf.Consequence.isin(intergenic_csqs).all():
        record.info['intergenic'] = True
        chrom = record.chrom
        if 'chr' not in chrom:
            chrom = 'chr' + chrom
        rec_bt = pbt.BedTool('{}\t{}\t{}\n'.\
                     format(chrom, record.pos, record.pos+1), 
                     from_string=True)
        try:
            dists = [int(f[-1]) for f in ras_bt.closest(rec_bt, d=True) if int(f[-1]) > -1]
            if len(dists) > 0:
                record.info['RAS_distance'] = int(np.nanmin(dists))
        except:
            msg = 'RAS_distance computation failed for record {}; skipping\n'
            stderr.write(msg.format(record.id))

    # Extract top ClinVar significance, if any reported
    if any(vdf.ClinVar_CLNSIG != ''):
        best_clinvar_idx = vdf.ClinVar_CLNSIG.map(clinvar_sig).sort_values().index.values[0]
        best_clinvar_label = clinvar_remap.get(vdf.ClinVar_CLNSIG[best_clinvar_idx])
        record.info['ClinVar'] = best_clinvar_label

    # Extract top COSMIC significance, if any reported
    if any(vdf.COSMIC_COSMIC_MUT_SIG != ''):
        best_cosmic = int(vdf.COSMIC_COSMIC_MUT_SIG.str.split('&', expand=True).\
                          max(axis=1).sort_values().head(1).values[0])
        record.info['COSMIC_mutation_tier'] = best_cosmic

    # Mve promoter annotations to INFO
    if any(vdf.promoters != ''):
        prom_genes = set(vdf.promoters.str.split(':', expand=True)[0])
        record.info['promoter'] = ','.join(prom_genes)

    # Move ABC enhancer annotations to INFO
    if any(vdf.ABC_enhancers != ''):
        record, vdf = _clean_abc(record, vdf)

    # Clean up SpliceAI (store in CSQ record matching gene)
    if any(vdf.SpliceAI_pred_SYMBOL != ''):
        record, vdf = _clean_spliceai(record, vdf, spliceai_cutoff)
        
    # Clean up values that have multiple scores per variant (MPC, FATHMM)
    def _multiMax(x, fx=np.nanmax):
        return fx([pd.to_numeric(k, errors='coerce') for k in x.split('&')])
    if any(vdf.MPC_score.str.contains('&')):
        vdf['MPC_score'] = vdf.MPC_score.map(_multiMax)
    if any(vdf.FATHMM_score.str.contains('&')):
        vdf['FATHMM_score'] = vdf.FATHMM_score.map(lambda x: _multiMax(x, fx=np.nanmin))

    # Clean up CCR and regional constraint (store in CSQ record matching gene)
    if any(vdf.CCR != ''):
        ccr = '&'.join(vdf.CCR[~vdf.CCR.duplicated()].values)
        ccr_dict = {x.split(':')[0] : float(x.split(':')[1]) for x in ccr.split('&')}
        vdf['CCR'] = vdf.SYMBOL.map(ccr_dict)
    if any(vdf.ExAC_regional_constraint != ''):
        rc = '&'.join(vdf.ExAC_regional_constraint[~vdf.ExAC_regional_constraint.duplicated()].values)
        rc_dict = {x.split(':')[1] : float(x.split(':')[3]) for x in rc.split('&')}
        vdf['ExAC_regional_constraint'] = vdf.SYMBOL.map(rc_dict)

    # Move GTEx eQTL information to INFO
    if any(vdf.GTEx != ''):
        record, vdf = _clean_gtex(record, vdf, tx_map)

    # Clean up Ensembl regulatory annotations
    if 'CELL_TYPE' in vdf.columns:
        record, vdf = _clean_regulatory(record, vdf, tx_map)

    # Overwrite CSQ INFO field with cleaned VEP info
    record.info.pop('CSQ')
    record.info['CSQ'] = tuple(['|'.join(vals.astype(str)) for i, vals in vdf.iterrows()])
    
    # Return cleaned record
    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input VEP-annotated .vcf')
    parser.add_argument('vcf_out', help='output .vcf [default: stdout]', default='stdout')
    parser.add_argument('--gtf', help='GTF for annotating CNVs', required=True)
    parser.add_argument('-t', '--transcript-info', required=True, help='.tsv ' + \
                        'mapping ENST:ENSG:symbol:length')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile('-', 'r')
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Parse mapping of VEP fields and transcript info
    vep_map = parse_vep_map(invcf)
    tx_map = load_tx_map(args.transcript_info)
    gtf_tabix = pysam.TabixFile(args.gtf)
    ras_bt = load_ras_bt(args.gtf)    
    
    # Open connection to output file
    out_header = reformat_header(invcf)
    if args.vcf_out in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile('-', 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=out_header)
    out_vep_map = parse_vep_map(outvcf)

    # Iterate over records in invcf, clean up each, and write to outvcf
    for record in invcf.fetch():
        if 'SVTYPE' in record.info.keys():
            record = cleanup_cnv(record, out_vep_map.values(), gtf_tabix)
        else:
            record = cleanup(record, vep_map, tx_map, ras_bt)
        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

