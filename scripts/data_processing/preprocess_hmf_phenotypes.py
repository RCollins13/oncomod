#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Preprocess HMF sample descriptives for relevant patients
"""


import argparse
import json
import pandas as pd
from datetime import date
from os.path import basename
from re import sub


def load_metadata(meta_in):
    """
    Load and clean HMF metadata.tsv
    """

    use_cols = 'sampleId setName hmfPatientId hmfSampleId tumorPurity gender ' + \
               'birthYear deathDate primaryTumorLocation primaryTumorSubLocation ' + \
               'primaryTumorType primaryTumorSubType biopsyDate ' + \
               'biopsySite biopsyLocation'

    meta = pd.read_csv(meta_in, sep="\t", usecols=use_cols.split())

    # Restrict to adenocarcinomas from CRC, lung, or pancreas
    keep_cloc = 'Colorectum Lung Pancreas Rectum Colon'.split()
    keep_ctype = 'Carcinoma Adenocarcinoma'.split()
    drop_stype = ['Small cell carcinoma', 'Neuroendocrine carcinoma',
                  'Squamous cell carcinoma']
    keep_idx = (meta.primaryTumorLocation.isin(keep_cloc)) & \
               ((meta.primaryTumorType.isin(keep_ctype)) | meta.primaryTumorType.isna()) & \
               (~meta.primaryTumorSubType.isin(drop_stype))
    meta = meta[keep_idx]

    # Clean or add simple variables as needed
    meta['gender'] = meta.gender.str.upper()
    meta['AJCC_STAGE'] = 4
    meta['APPROX_STAGE'] = 4
    meta['CAUSE_OF_DEATH'] = pd.NA
    meta['LAST_ALIVE_DATE'] = pd.NA
    meta['BIOPSY_SITE_TYPE'] = pd.NA
    for i in range(10):
        meta['PC' + str(i+1)] = pd.NA
    meta['POPULATION'] = pd.NA
    meta['NORMAL_ID'] = meta.sampleId.apply(lambda x: sub('T[0-9]*$', 'R', x))

    # Convert cancer labels to study abbreviations
    ctype_map = {'Colorectum' : 'CRAD',
                 'Rectum' : 'CRAD',
                 'Colon' : 'CRAD',
                 'Lung' : 'LUAD',
                 'Pancreas' : 'PDAC'}
    meta['CANCER_TYPE'] = meta.primaryTumorLocation.map(ctype_map)

    # Generate other missing labels
    def _build_diagnosis(vals):
        if not pd.isna(vals.primaryTumorSubLocation):
            p1 = vals.primaryTumorSubLocation
        elif not pd.isna(vals.primaryTumorLocation):
            p1 = vals.primaryTumorLocation
        else:
            p1 = None
        if not pd.isna(vals.primaryTumorSubType):
            p2 = vals.primaryTumorSubType.lower()
        elif not pd.isna(vals.primaryTumorType):
            p2 = vals.primaryTumorType.lower()
        else:
            p2 = None
        return ' '.join([str(v) for v in [p1, p2] if v is not None])
    meta['PRIMARY_CANCER_DIAGNOSIS'] = meta.apply(_build_diagnosis, axis=1)
    meta['IS_ALIVE'] = meta.deathDate.isna().astype(int)

    # Manipulate all date-based fields
    meta['DIAGNOSIS_YEAR'] = meta.biopsyDate.map(lambda x: str(x).split('-')[0]).\
                                             astype('float').astype('Int64')
    meta['AGE_AT_DIAGNOSIS'] = meta['DIAGNOSIS_YEAR'] - meta.birthYear.astype("Int64")
    def _days_survived(vals):
        bd = vals.biopsyDate
        dd = vals.deathDate
        if any(pd.isna([bd, dd])):
            return pd.NA
        else:
            bd = date(*[int(v) for v in bd.split('-')])
            dd = date(*[int(v) for v in dd.split('-')])
            return (dd - bd).days
    meta['DAYS_SURVIVED'] = meta.apply(_days_survived, axis=1)

    return meta


def polish_metadata(meta, purple_in):
    """
    Polish missing sample metadata using Purple somatic QC information
    Also excludes MSI tumors
    """

    pdf = pd.read_csv(purple_in, sep='\t').rename(columns={'#sample' : 'sampleId'})

    # Subset to microsatellite stable tumors
    mss_ids = pdf.loc[pdf.msStatus == 'MSS', 'sampleId']
    meta = meta[meta.sampleId.isin(mss_ids)]

    # Fill missing sex labels based on Purple estimated sex
    nosex = meta.gender.isna()
    psex = pdf.set_index('sampleId').gender.to_dict()
    meta.loc[nosex, 'gender'] = meta.loc[nosex, 'sampleId'].map(psex)

    return meta


def select_patients(meta, vcf_paths=None, json_in=None):
    """
    Subset to unique patients with matching germline data (if optioned)
    """

    if vcf_paths is not None and json_in is not None:

        # Infer sample IDs based on file paths
        with open(vcf_paths) as fin:
            vcfs = [basename(l.rstrip()) for l in fin.readlines()]
        sids = [v.split('.')[0] for v in vcfs]

        # Get mapping of tumor IDs to germline VCF prefixes from json
        vcf_prefix_map = {}
        with open(json_in) as fin:
            mfst = json.load(fin)
            for sdict in mfst['data']:
                sid = sdict['sampleId']
                vcf = sdict['germline']['variants']['url']
                vcf_prefix = basename(vcf).split('.')[0]
                vcf_prefix_map[sid] = vcf_prefix

        # Only keep rows with matching germline VCFs
        keepers = meta.sampleId.map(vcf_prefix_map).isin(sids)
        meta = meta[keepers]

    # Deduplicate based on unique patient IDs
    meta.drop_duplicates(subset='hmfPatientId', inplace=True)

    return meta


def clean_output_df(meta):
    """
    Clean up main patient data frame before writing to output file
    """

    new_col_names = {'NORMAL_ID' : 'SAMPLE_ID',
                     'gender' : 'SEX',
                     'biopsyDate' : 'DIAGNOSIS_DATE',
                     'biopsySite' : 'BIOPSY_SITE',
                     'tumorPurity' : 'TUMOR_PURITY'}
    drop_cols = 'sampleId setName hmfPatientId hmfSampleId birthYear ' + \
                'deathDate primaryTumorLocation primaryTumorSubLocation ' + \
                'primaryTumorType primaryTumorSubType'

    meta.rename(columns=new_col_names, inplace=True)
    meta.drop(drop_cols.split(), axis=1, inplace=True)
    col_order = 'SAMPLE_ID SEX POPULATION CANCER_TYPE PRIMARY_CANCER_DIAGNOSIS ' + \
                'AJCC_STAGE APPROX_STAGE DIAGNOSIS_DATE AGE_AT_DIAGNOSIS ' + \
                'DIAGNOSIS_YEAR IS_ALIVE LAST_ALIVE_DATE DAYS_SURVIVED ' + \
                'CAUSE_OF_DEATH PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 ' + \
                'BIOPSY_SITE BIOPSY_SITE_TYPE TUMOR_PURITY'
    sort_cols = 'CANCER_TYPE AGE_AT_DIAGNOSIS POPULATION'

    return meta[col_order.split()].sort_values(sort_cols.split())


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metadata', help='HMF metadata.tsv', required=True)
    parser.add_argument('--purple-qc', help='.tsv of Purple somatic QC')
    parser.add_argument('--germline-vcfs-list', help='List of germline VCFs')
    parser.add_argument('--hmf-json', help='HMF manifest.json')
    parser.add_argument('--out-prefix', help='Prefix for output files', required=True)
    args = parser.parse_args()

    # Load & clean metadata
    meta = load_metadata(args.metadata)

    # Polish metadata using Purple QC
    if args.purple_qc is not None:
        meta = polish_metadata(meta, args.purple_qc)

    # Restrict to tumors from patients with germline data available
    meta = select_patients(meta, args.germline_vcfs_list, args.hmf_json)

    # Add ancestry information
    # TODO: need to implement this once we have PCs

    # Clean up columns for final outputs
    meta = clean_output_df(meta)

    import pdb; pdb.set_trace()

    # Write out sample IDs per cancer type
    for cancer in 'PDAC CRAD LUAD'.split():
        with open(args.out_prefix + cancer + '.samples.list', 'w') as fout:
            ids = set(meta.DONOR_ID.astype(str)[meta.CANCER_TYPE == cancer].to_list())
            for sample in ids:
                fout.write(sample + '\n')
    all_cancers_out = args.out_prefix + 'ALL' + '.samples.with_cancer_labels.tsv'
    meta.loc[:, 'DONOR_ID CANCER_TYPE'.split()].\
            to_csv(all_cancers_out, sep='\t', index=False, header=False)
    all_samples_out = args.out_prefix + 'ALL' + '.samples.list'
    meta.DONOR_ID.to_csv(all_samples_out, sep='\t', index=False, header=False)

    # Write out full table of patient metadata for all cancer types
    meta.rename(columns = {'DONOR_ID' : '#DONOR_ID'}, inplace=True)
    meta.to_csv(args.out_prefix + 'ALL.sample_metadata.tsv.gz', 
                   index=False, sep='\t', na_rep='.')


if __name__ == '__main__':
    main()
