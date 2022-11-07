#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Subset DFCI-Profile somatic data to relevant tumors
"""


import argparse
import numpy as np
import pandas as pd
import pysam
from sys import stdout, stderr


def load_id_map(id_map_in):
    """
    Load mappings between exome, array, and donor IDs
    """

    id_df = pd.read_csv(id_map_in, sep='\t').rename(columns={'#DONOR_ID' : 'DONOR_ID'})

    return {'exome_to_array' : id_df.set_index(id_df.WES_BAM_ID).ARRAY_TYPED_ID.to_dict(),
            'array_to_exome' : id_df.set_index(id_df.ARRAY_TYPED_ID).WES_BAM_ID.to_dict(),
            'exome_to_donor' : id_df.set_index(id_df.WES_BAM_ID).DONOR_ID.to_dict(),
            'array_to_donor' : id_df.set_index(id_df.ARRAY_TYPED_ID).DONOR_ID.to_dict()}


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


def determine_next_record(ex_rec, gt_rec, imp_rec):
    """
    Determine which record(s) to process & write to file next
    """

    cur_pos = {}
    for i, rec in enumerate([ex_rec, gt_rec, imp_rec]):
        if rec is None:
            cur_pos[i] = np.nan
        else:
            cur_pos[i] = rec.pos
    cur_pos = pd.Series(cur_pos)

    next_idxs = np.where(cur_pos == np.nanmin(cur_pos))[0].tolist()

    return [k for i, k in enumerate('exome array-typed array-imputed'.split()) if i in next_idxs]


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


def resolve_array_discrepancies(current_records, current_ids, in_vcfs, next_techs, 
                                ref_fa, verbose=False):
    """
    Main wrapper for handling discrepancies between typed and imputed array records
    """

    # First, assess GT correlation between typed and imputed samples
    if current_ids['array-typed'] != current_ids['array-imputed']:
        gt_cor = calc_gt_correlation(current_records['array-typed'],
                                     current_records['array-imputed'])

        # If genotypes are inversely correlated between typed & imputed,
        # assume that array probe incorrectly targeted ref. allele as "alt.".
        # In this case, can use typed GTs but need to reverse ref/alt
        # and GTs for all samples accordingly
        if gt_cor <= -0.5:
            current_records['array-typed'], current_ids['array-typed'] = \
                invert_genotypes(current_records['array-typed'], gt_cor, verbose)

    # Next, check for strand flips in the array typed record
    # To be conservative, this must require imputed ref allele to match
    # reference fasta allele
    if current_ids['array-typed'] != current_ids['array-imputed']:
        ref_allele = ref_fa.fetch(current_records['array-imputed'].chrom, 
                                  current_records['array-imputed'].pos,
                                  current_records['array-imputed'].pos+1)
        if current_records['array-imputed'].alleles[0] == ref_allele:
            current_records['array-typed'], current_ids['array-typed'] = \
                correct_strand(current_records['array-typed'], ref_fa, verbose)


    # If neither GT inversion or strand flip corrects discrepancies, keep the
    # record that matches the canonical ref allele
    if current_ids['array-typed'] != current_ids['array-imputed']:
        ref_allele = ref_fa.fetch(current_records['array-imputed'].chrom, 
                                  current_records['array-imputed'].pos,
                                  current_records['array-imputed'].pos+1)
        skip_techs = []
        if current_records['array-typed'].alleles[0] != ref_allele:
            skip_techs.append('array-typed')
        if current_records['array-imputed'].alleles[0] != ref_allele:
            skip_techs.append('array-imputed')
        if verbose:
            if len(skip_techs) == 1:
                keep_tech = list(set('array-typed array-imputed'.split()).difference(set(skip_techs)))[0]
                fmt = '  ** Could not resolve discrepancies between typed and imputed ' + \
                      ' array records at {:,}. Keeping {} as it matches the reference\n'
                stderr.write(fmt.format(current_records[keep_tech].pos, keep_tech))
            else:
                fmt = '  ** Could not resolve discrepancies between typed and imputed ' + \
                      ' array records at {:,}. Skipping both as neither match the reference\n'
                stderr.write(fmt.format(current_records['array-typed'].pos))
                dbg = True
        for skip_tech in skip_techs:
            next_techs = [t for t in next_techs if t != skip_tech]
            new_rec, new_id = next_record(in_vcfs[skip_tech])
            current_records[skip_tech] = new_rec
            current_ids[skip_tech] = new_id

    return current_records, current_ids, in_vcfs, next_techs


def map_genotypes(old_rec, new_rec, id_mappings, rec_fmt):
    """
    Extracts GTs in old_rec and maps them onto samples in new_rec
    Assumes id_mappings is dict keyed on samples in old_rec
    Annotates AC, AN, and AF while mapping
    """

    AN, AC, AF = [0] * 3

    for sample in old_rec.samples:

        # Get old genotype
        # Note: DeepVariant (used on exomes) does not call ref GTs
        # We have already pre-filtered for well-captured exome intervals,
        # so any no-call GT from exomes should be treated as 0/0
        GT = old_rec.samples[sample]['GT']
        if GT == (None, None):
            if rec_fmt == 'exome':
                GT = (0, 0)
        else:
            if None in GT:
                GT = (None, [g for g in GT if g is not None][0])
            else:
                GT = tuple(sorted(GT))

        # Assign genotype to new record
        new_rec.samples[id_mappings[sample]]['GT'] = GT

        # Update AC/AN
        AC += len([a for a in GT if a is not None and a > 0])
        AN += len([a for a in GT if a is not None])

    # Annotate AC/AN/AF in record
    new_rec.info['AC'] = AC
    new_rec.info['AN'] = AN
    if AN > 0:
        new_rec.info['AF'] = AC / AN
    
    return new_rec


def format_record(old_rec, id_map, out_vcf, rec_fmt):
    """
    Reformat an existing variant record and write to out_vcf
    """

    # Copy basic information from old record to new record
    new_id = '_'.join([rec_fmt.replace('-', '_'), old_rec.chrom, 
                       str(old_rec.pos), *old_rec.alleles])
    new_rec = out_vcf.new_record(contig=old_rec.chrom, start=old_rec.pos, 
                                 alleles=old_rec.alleles, id=new_id)
    new_rec.info['SOURCE'] = rec_fmt.replace('-', '_')
    new_rec.info.pop('END')

    # Map genotypes
    if rec_fmt == 'exome':
        new_rec = map_genotypes(old_rec, new_rec, id_map['exome_to_donor'], rec_fmt)
    else:
        new_rec = map_genotypes(old_rec, new_rec, id_map['array_to_donor'], rec_fmt)

    out_vcf.write(new_rec)


def report_current_coords(current_records):
    """
    Print position of current records
    """

    fmt = 'Current positions (EX|GT|IMP): {}\n'
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
    parser.add_argument('--sample-id-map', help='.tsv mapping between various IDs', required=True)
    parser.add_argument('--exome-vcf', help='Exome .vcf', required=True)
    parser.add_argument('--array-typed-vcf', help='Array genotyped .vcf', required=True)
    parser.add_argument('--array-imputed-vcf', help='Array imputed .vcf', required=True)
    parser.add_argument('--ref-fasta', help='Reference .fasta', required=True)
    parser.add_argument('--header', help='Header to use for output .vcf', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to merged ' +
                        '.vcf [default: stdout]')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable ' +
                        'verbose logging [default: no logging]')
    args = parser.parse_args()

    # Load sample ID map
    id_map = load_id_map(args.sample_id_map)

    # Open connection to reference fa
    ref_fa = pysam.FastaFile(args.ref_fasta)

    # Load header and add all donor IDs
    header = pysam.VariantFile(args.header).header
    for sample in id_map['exome_to_donor'].values():
        header.add_sample(sample)

    # Open connections to all three input VCFs
    in_vcfs = {'exome' : pysam.VariantFile(args.exome_vcf),
               'array-typed' : pysam.VariantFile(args.array_typed_vcf),
               'array-imputed' : pysam.VariantFile(args.array_imputed_vcf)}

    # Open connection to output VCF
    fout = args.outfile
    if fout in '- stdout /dev/stdout'.split():
        fout = stdout
    out_vcf = pysam.VariantFile(args.outfile, 'w', header=header)

    # Walk through all three VCFs based on lowest coordinate
    variants_seen = set()
    ex_rec, ex_id = next_record(in_vcfs['exome'])
    gt_rec, gt_id = next_record(in_vcfs['array-typed'])
    imp_rec, imp_id = next_record(in_vcfs['array-imputed'])
    current_records = {'exome' : ex_rec, 'array-typed' : gt_rec, 'array-imputed' : imp_rec}
    current_ids = {'exome' : ex_id, 'array-typed' : gt_id, 'array-imputed' : imp_id}

    # Keep processing records until all VCFs are exhausted
    while not all(r is None for r in current_records.values()):

        if args.verbose:
            report_current_coords(current_records)

        # Determine which record(s) to process next
        next_techs = determine_next_record(*current_records.values())

        # Resolve discrepancies between genotyped & imputed data
        if 'array-typed' in next_techs \
        and 'array-imputed' in next_techs \
        and current_ids['array-typed'] != current_ids['array-imputed']:
            current_records, current_ids, in_vcfs, next_techs = \
                resolve_array_discrepancies(current_records, current_ids, in_vcfs, 
                                            next_techs, ref_fa, args.verbose)

        # Process current file pointer(s) that have the left-most position
        for tech in next_techs:

            # Write record to VCF if it hasn't already been written by a 
            # higher-priority technology
            if current_ids[tech] not in variants_seen:
                format_record(current_records[tech], id_map, out_vcf, tech)
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

