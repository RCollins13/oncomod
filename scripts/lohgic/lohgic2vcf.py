#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert final germline predictions from PROFILE to dense VCF format
"""


import argparse
import pandas as pd
import pybedtools as pbt
import pysam
from datetime import datetime
from math import ceil, log10
from numpy import nanmedian
from sys import stdout


missing_flag_base = 'NOT_CAPTURED_IN_ONCOPANEL_v{}'


def report_progress(n, k, step, then, begin):
    """
    Report progress
    """

    now = datetime.now()
    tdelta = (now - then).seconds
    ttotal = (now - begin).seconds
    spv = ttotal / k
    then = now
    tremain = ceil((spv * (n - k)) / 60)
    msg = 'lohgic2vcf.py: Processed {:,} records in ~{:,} seconds'
    print(msg.format(step, tdelta))
    msg = '               * Rolling average: {:.2f} seconds / record'
    print(msg.format(spv))
    msg = '               * Currently {:.2f}% complete ({:,} records finished of {:,} total)'
    print(msg.format(100 * k / n, k, n))
    msg = '               * Approximately {} minutes remaining\n'
    print(msg.format(tremain))

    return then


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-L', '--lohgic', help='path to .tsv output by ' +
                        'germline_classifier.R', required=True)
    parser.add_argument('-H', '--vcf-header', help='header for output .vcf',
                        required=True)
    parser.add_argument('-S', '--samples-list', help='list of sample IDs', 
                        required=True)
    parser.add_argument('-P', '--pbp-map', help='.tsv mapping MRNs and BL IDs ' +
                        'to PBPs', required=True)
    parser.add_argument('-i', '--capture-intervals', help='.bed with intervals ' +
                        'captured per oncopanel version (fourth column). ' +
                        'If provided, genotypes from samples not captured at a ' +
                        'position will be replaced by no-calls')
    parser.add_argument('-m', '--seq-metrics', help='.tsv mapping seq depth to biopsy ID. ' +
                        'If provided, used to estimate reference genotype quality.')
    parser.add_argument('-o', '--vcf-out', help='path to output .vcf ' +
                        '[default: stdout]', default='stdout')
    parser.add_argument('-p', '--report-progress', action='store_true',
                        help='Report progress')
    args = parser.parse_args()

    # Load labeled LOHGIC data as pd.DataFrame
    mdf = pd.read_csv(args.lohgic, sep='\t', low_memory=False)

    # Retain only germline predictions
    mdf = mdf[mdf.pred_label == 'germline']

    # Sort variants by chrom, pos, ref, len(alt), alt
    mdf.loc[:, 'alt_len'] = mdf.ALT_ALLELE.fillna('N').map(len) - mdf.REF_ALLELE.fillna('N').map(len)
    mdf.sort_values('CHROMOSOME POSITION REF_ALLELE alt_len ALT_ALLELE'.split(), inplace=True)

    # Compute GQs for all variants
    mdf.loc[:, 'GQ'] = -10 * (1 - mdf.prob_germline).map(log10)

    # Open connection to output VCF
    out_header = pysam.VariantFile(args.vcf_header).header
    if args.vcf_out in '- stdout /dev/stdout':
        outpath = stdout
    else:
        outpath = args.vcf_out
    out_vcf = pysam.VariantFile(outpath, mode='w', header=out_header)

    # Load list of all samples (expected to be PBPs)
    all_samples = set()
    with open(args.samples_list) as fin:
        for sid in fin.readlines():
            sid = sid.rstrip()
            all_samples.add(sid)
            out_vcf.header.samples.add(sid)

    # Load PBP map and convert all BL IDs in mdf to PBPs
    pbpdf = pd.read_csv(args.pbp_map, sep='\t')
    pbpdf.index = pbpdf.SAMPLE_ACCESSION_NBR
    pbp_map = pbpdf.PBP.to_dict()
    mdf.loc[:, 'PBP'] = mdf.SAMPLE_ACCESSION_NBR.map(pbp_map)

    # Load captured interval data, if optioned
    if args.capture_intervals is not None:
        capture_intervals = pbt.BedTool(args.capture_intervals)

        # Also build map of sample ID to panel version
        pbpdf.index = pbpdf.PBP
        panel_map = pbpdf.PANEL_VERSION.to_dict()
        samples_by_version = {}
        for pbp, v in panel_map.items():
            if v not in samples_by_version.keys():
                samples_by_version[v] = set([pbp])
            else:
                samples_by_version[v].add(pbp)

        # Also update header with INFO for uncaptured flags
        info_base = '##INFO=<ID={},Number=0,Type=Flag,Description="This site ' + \
                    'was not captured in OncoPanel version {}">'
        for v in set(panel_map.values()):
            out_vcf.header.add_line(info_base.format(missing_flag_base.format(v), v))

    # Load sample depth and convert to approx ref GQs
    if args.seq_metrics is not None:
        seqdf = pd.read_csv(args.seq_metrics, sep='\t')
        seqdf.index = seqdf.SAMPLE_ACCESSION_NBR.map(pbp_map)
        seqdf.loc[:, 'GQ'] = -10 * (1 / seqdf.MEAN_SAMPLE_COVERAGE).map(log10)
        ref_gq = seqdf.GQ.to_dict()
    else:
        ref_gq = {}

    # Prepare for progress logging, if optioned
    gb = mdf.groupby(by='CHROMOSOME POSITION REF_ALLELE ALT_ALLELE'.split())
    if args.report_progress:
        n = len(gb)
        k = 0
        step = 100
        msg = 'lohgic2vcf.py: Beginning to convert {:,} unique variants ' + \
              'in input .tsv to VCF. Please be patient...\n'
        print(msg.format(n))
        begin = datetime.now()
        then = begin

    # Iterate over mdf grouped by chrom, pos, ref, alt
    for gbi, subdf in gb:

        # Set basic record info
        record = out_vcf.new_record()
        chrom, pos, ref, alt = gbi
        record.chrom = str(chrom)
        record.pos = int(pos)
        record.id = '_'.join([str(k) for k in gbi])
        record.alleles = (ref, alt)
        record.filter.add('PASS')

        # Check for capture intervals, if optioned
        if args.capture_intervals is not None:
            rec_bt = pbt.BedTool('{}\t{}\t{}\n'.format(*gbi[:2], gbi[1] + 1), 
                                 from_string=True)
            captured_in = set()
            for hit in capture_intervals.intersect(rec_bt, wa=True):
                captured_in.add(float(hit[3]))
            for v in '1 2 3 3.1'.split():
                if float(v) not in captured_in:
                    record.info[missing_flag_base.format(float(v))] = True
        
        # Get list of non-reference PBPs and their GQs
        subdf.index = subdf.PBP
        nonref_gq = subdf.GQ.to_dict()
        if len(nonref_gq) == 0:
            if args.report_progress:
                k += 1
                if k % step == 0:
                    then = report_progress(n, k, step, then, begin)
            continue

        # Set QUAL as the median non-ref GQ
        record.qual = round(nanmedian(list(nonref_gq.values())), 2)

        # Iterate over all PBPs and write genotypes
        AC, AN = 0, 0
        for pbp in all_samples:

            # If sample is non-reference, add to GT
            if pbp in nonref_gq.keys():
                record.samples[pbp]['GT'] = (0, 1)
                record.samples[pbp]['GQ'] = round(nonref_gq.get(pbp, None), 2)
                AC += 1
                AN += 1

            else:
                # If --capture-intervals is not provided, assume 0/0 GT
                if args.capture_intervals is None:
                    record.samples[pbp]['GT'] = (0, 0)
                    record.samples[pbp]['GQ'] = round(ref_gq.get(pbp, None), 2)
                    AN += 1
                    continue

                # Otherwise, check if this position was captured in this sample
                version = float(panel_map.get(pbp, None))
                if version in captured_in:
                    record.samples[pbp]['GT'] = (0, 0)
                    record.samples[pbp]['GQ'] = round(ref_gq.get(pbp, None), 2)
                    AN += 1
                else:
                    record.samples[pbp]['GT'] = (None, None)
                    record.samples[pbp]['GQ'] = None

        # Update allele frequency info
        record.info['AN'] = AN
        record.info['AC'] = AC
        if AN > 0:
            AF = AC / AN
        else:
            AF = None
        record.info['AF'] = AF

        # Write record to output VCF
        if AC > 0:
            out_vcf.write(record)

        # Report progress, if optioned
        if args.report_progress:
            k += 1
            if k % step == 0:
                then = report_progress(n, k, step, then, begin)


    # Close connection to output VCF to clear buffer
    out_vcf.close()


if __name__ == '__main__':
    main()

