#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Define a subset of "control" samples for germline-somatic association based on
their tumor alterations
"""


import argparse
import json
import pandas as pd
import pybedtools as pbt
import pysam
from sys import stdout, path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from vep_utils import parse_vep_map, vep2df


def check_rules(vdf, rules):
    """
    Checks if a variant meets _all_ of the criteria specified in rules
    `rules` should be a dict mapping VEP keys to desired values

    Default return is True unless any rules are failed
    """

    hit_idxs = pd.Series([True] * len(vdf), index=vdf.index)

    for key, value in rules.items():

        if key == 'Amino_acids':
            hit_idxs = vdf[key].str.split('/').map(lambda x: value in x)

        elif key in vdf.columns:
            hit_idxs = vdf.loc[hit_idxs, key].astype(str) == str(value)

        else:
            print('\nNot sure how to interpret key {}\n'.format(key))
            import pdb; pdb.set_trace()

        if not hit_idxs.any():
            break

    return hit_idxs.any()


def check_rules_list(vdf, rules_list):
    """
    Check whether a VEP DataFrame meets any of the criteria sets in rules_list
    rules_list should be a list of dicts mapping VEP key to desired value
    """

    # Check compliance for each set of rules in rules_list
    hits = [check_rules(vdf, rset) for rset in rules_list]
    
    return any(hits)


def check_all_criteria(record, criteria, vep_map):
    """
    Check if a record meets any criteria
    """

    hit = False

    vdf = vep2df(record, vep_map)

    for gene, rules in criteria.items():

        gene_rows = (vdf.SYMBOL == gene)
        if not gene_rows.any():
            continue

        hit = check_rules_list(vdf[gene_rows], rules)

        if hit:
            break

    return hit


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='input vcf.')
    parser.add_argument('--criteria', required=True, help='json of exclusion ' + 
                        'criteria; samples meeting any of these criteria will ' +
                        'not be included in --outfile. Criteria format should be ' +
                        'JSON object keyed by gene symbol, where each element in ' +
                        'the object is a list of criteria referring to a specific ' +
                        'mutation/alteration/codon/etc')
    parser.add_argument('--regions', help='Subset search of --vcf to regions ' + 
                        'specified in this BED file.')
    parser.add_argument('--eligible-samples', help='list of eligible samples. If ' +
                        'not provided, all samples in --vcf will be used.')
    parser.add_argument('--exclude-samples', help='list of samples to exclude ' +
                        'from all evaluations')
    parser.add_argument('-o', '--outfile', help='output .txt [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Load criteria
    with open(args.criteria) as fin:
        criteria = json.load(fin)

    # Open connection to input vcf
    vcf = pysam.VariantFile(args.vcf)

    # Get mapping of VEP CSQ fields
    vep_map = parse_vep_map(vcf)

    # Build list of all eligible samples
    if args.eligible_samples is not None:
        with open(args.eligible_samples) as fin:
            elig_samples = set([l.rstrip() for l in fin.readlines()])
    else:
        elig_samples = set(vcf.header.samples)
    if args.exclude_samples is not None:
        with open(args.exclude_samples) as fin:
            excl_samples = set([l.rstrip() for l in fin.readlines()])
            elig_samples = elig_samples.difference(excl_samples)

    # Load regions to fetch
    if args.regions is not None:
        regions = pbt.BedTool(args.regions)
    else:
        regions = []

    # Iterate over records in VCF
    for region in regions:

        if len(region) > 0:
            fetcher = vcf.fetch(region.chrom, region.start, region.end)
        else:
            fetcher = vcf.fetch()
        
        for record in fetcher:

            # Check if record meets any control exclusion criteria
            hit = check_all_criteria(record, criteria, vep_map)

            # If record meets criteria, remove all non-ref samples from list of eligible samples
            if hit:
                for sid, svals in record.samples.items():
                    gt = [a for a in svals['GT'] if a is not None]
                    if sum(gt) > 0:
                        elig_samples.discard(sid)

    # Write results to output tsv
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    for sid in sorted(list(elig_samples)):
        outfile.write(sid + '\n')
    outfile.close()


if __name__ == '__main__':
    main()

