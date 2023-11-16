#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Merge imputed SNPs and predicted coding variants for DFCI-PROFILE cohort
"""


import argparse
import numpy as np
import pandas as pd
import pysam
from sys import stdout, stderr


def next_record(vcf):
    """
    Helper to advance VCF by one record and return a unique, hashable variant ID string
    """

    try:
        rec = vcf.__next__()
        vid = '{}_{}_{}_{}'.format(rec.chrom, rec.pos, *rec.alleles)
    except:
        rec, vid = None, None

    return rec, vid


def determine_next_record(imp_rec, ex_rec):
    """
    Determine which record(s) to process & write to file next
    """

    cur_pos = {}
    for i, rec in enumerate([imp_rec, ex_rec]):
        if rec is None:
            cur_pos[i] = np.nan
        else:
            cur_pos[i] = rec.pos
    cur_pos = pd.Series(cur_pos)

    next_idxs = np.where(cur_pos == np.nanmin(cur_pos))[0].tolist()

    return [k for i, k in enumerate('imputed lohgic'.split()) if i in next_idxs]


def correct_strand(record, fasta, verbose=False):
    """
    Corrects strand flips in genotyped array data
    """

    # Check if record is a candidate for a strand flip
    old_alleles = record.alleles
    if len(record.alleles) == 2:
        ref_allele = fasta.fetch(record.chrom, record.pos, record.pos+1)
        if ref_allele == record.alleles[1]:
            record.alleles = record.alleles[::-1]

            if verbose:
                fmt = '  ** Detected strand flip at {:,} (typed as {} / {} but should be {} / {})\n'
                stderr.write(fmt.format(record.pos, *old_alleles, *record.alleles))

    vid = '{}_{}_{}_{}'.format(record.chrom, record.pos, *record.alleles)

    return record, vid


def calc_gt_correlation(rec1, rec2):
    """
    Compute simple correlation coefficient between two sets of genotypes
    """

    samp1 = set(s for s in rec1.samples if rec1.samples[s]['GT'] is not (None, None))
    samp2 = set(s for s in rec2.samples if rec2.samples[s]['GT'] is not (None, None))
    samples = samp1.intersection(samp2)

    def _get_ad(gt):
        return sum([a for a in gt if a is not None])

    ad1 = np.array([_get_ad(rec1.samples[s]['GT']) for s in samples])
    ad2 = np.array([_get_ad(rec2.samples[s]['GT']) for s in samples])

    return np.corrcoef(ad1, ad2)[0, 1]


def invert_genotypes(record, gt_cor, verbose=False):
    """
    Inverts ref/alt and all genotypes for a single record
    """

    nt_map = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    new_alleles = tuple([nt_map[a] for a in record.alleles])


    if verbose:
        fmt = '  ** Detected inverted record at {:,} (typed as {} / {} but should be {} / {}; r={:.2f})\n'
        stderr.write(fmt.format(record.pos, *record.alleles, *new_alleles, gt_cor))

    record.alleles = new_alleles

    for sample in record.samples.keys():
        old_gt = record.samples[sample]['GT']
        if old_gt != (None, None):
            record.samples[sample]['GT'] = tuple([abs(a - 1) for a in old_gt])

    vid = '{}_{}_{}_{}'.format(record.chrom, record.pos, *record.alleles)

    return record, vid


def map_genotypes(old_rec, new_rec, rec_fmt):
    """
    Extracts GTs in old_rec and maps them onto samples in new_rec
    Assumes id_mappings is dict keyed on samples in old_rec
    Annotates AC, AN, and AF while mapping
    """

    AN, AC, AF = [0] * 3

    old_samples = old_rec.samples.keys()

    for sample in new_rec.samples.keys():

        # Get old genotype
        if sample in old_samples:
            GT = old_rec.samples[sample]['GT']
            if GT != (None, None, ):
                GT = tuple(sorted(GT))
        else:
            GT = (None, None, )

        # Assign genotype to new record
        new_rec.samples[sample]['GT'] = GT

        # Update AC/AN
        AC += len([a for a in GT if a is not None and a > 0])
        AN += len([a for a in GT if a is not None])

    # Annotate AC/AN/AF in record
    new_rec.info['AC'] = AC
    new_rec.info['AN'] = AN
    if AN > 0:
        new_rec.info['AF'] = AC / AN
    
    return new_rec


def format_record(old_rec, out_vcf, rec_fmt):
    """
    Reformat an existing variant record and write to out_vcf
    """

    # Copy basic information from old record to new record
    new_id = '_'.join([rec_fmt.replace('-', '_'), old_rec.chrom, 
                       str(old_rec.pos), *old_rec.alleles])
    new_rec = out_vcf.new_record(contig=old_rec.chrom, start=old_rec.pos-1, 
                                 alleles=old_rec.alleles, id=new_id)
    # Note: subtracting one from old_rec.pos is necessary because pysam automatically
    # increments pos by one (assuming this is related to 0- vs 1- indexing but confirmed 
    # manually this is necessary to match reference nucleotides)
    new_rec.info['SOURCE'] = rec_fmt.replace('-', '_')
    new_rec.info.pop('END')

    # Populate genotypes
    new_rec = map_genotypes(old_rec, new_rec, rec_fmt)

    out_vcf.write(new_rec)


def report_current_coords(current_records):
    """
    Print position of current records
    """

    fmt = 'Current positions (IMPUTED|LOHGIC): {}\n'
    positions = []
    for rec in current_records.values():
        if rec is None:
            positions.append('-')
        else:
            positions.append('{:,}'.format(rec.pos))
    stderr.write(fmt.format(' | '.join(positions)))


def report_variant_processed(current_records, tech, outcome):
    """
    Report the results of processing a single record
    """

    fmt = '  - Record from {} {}: {:,} ({} / {})\n'
    rec = current_records[tech]
    stderr.write(fmt.format(tech, outcome, rec.pos, *rec.alleles))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--lohgic-vcf', help='LOHGIC germline .vcf', required=True)
    parser.add_argument('--imputed-vcf', help='Imputed .vcf', required=True)
    parser.add_argument('--header', help='Header to use for output .vcf', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to merged ' +
                        '.vcf [default: stdout]')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable ' +
                        'verbose logging [default: no logging]')
    args = parser.parse_args()

    # Open connections to all three input VCFs
    in_vcfs = {'imputed' : pysam.VariantFile(args.imputed_vcf),
               'lohgic' : pysam.VariantFile(args.lohgic_vcf)}

    # Gather list of all samples present in either VCF
    all_samples = set()
    for f in in_vcfs.values():
        all_samples.update(set([s for s in f.header.samples]))
    all_samples = sorted(list(all_samples))

    # Load header and add all donor IDs
    header = pysam.VariantFile(args.header).header
    for sample in all_samples:
        header.add_sample(sample)

    # Open connection to output VCF
    fout = args.outfile
    if fout in '- stdout /dev/stdout'.split():
        fout = stdout
    out_vcf = pysam.VariantFile(args.outfile, 'w', header=header)

    # Walk through both VCFs based on lowest coordinate
    variants_seen = set()
    imp_rec, imp_id = next_record(in_vcfs['imputed'])
    ex_rec, ex_id = next_record(in_vcfs['lohgic'])
    current_records = {'imputed' : imp_rec, 'lohgic' : ex_rec}
    current_ids = {'imputed' : imp_id, 'lohgic' : ex_id}

    # Keep processing records until all VCFs are exhausted
    while not all(r is None for r in current_records.values()):

        if args.verbose:
            report_current_coords(current_records)

        # Determine which record(s) to process next
        next_techs = determine_next_record(*current_records.values())

        # Process current file pointer(s) that have the left-most position
        for tech in next_techs:

            # Write record to VCF if it hasn't already been written by a 
            # higher-priority technology
            if current_ids[tech] not in variants_seen:
                format_record(current_records[tech], out_vcf, tech)
                variants_seen.add(current_ids[tech])
                if args.verbose:
                    report_variant_processed(current_records, tech, 'written')
            else:
                if args.verbose:
                    report_variant_processed(current_records, tech, 'skipped')

            # Advance this technology to the next record
            new_rec, new_id = next_record(in_vcfs[tech])
            current_records[tech] = new_rec
            current_ids[tech] = new_id

    # Close outfile to clear buffer
    out_vcf.close()


if __name__ == '__main__':
    main()

