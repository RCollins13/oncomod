#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Apply user-specified set(s) of rule(s) to a VCF to define variant sets for burden testing
"""


import argparse
import json
import numpy as np
import operator
import os
import pandas as pd
import pysam
import re
from sys import stdin, stdout, path
path.insert(0, os.path.join(path[0], '..', '..', 'utils'))
from vep_utils import parse_vep_map, vep2df


ops = {'==' : operator.eq,
       '!=' : operator.ne,
       '<' : operator.lt,
       '<=' : operator.le,
       '>' : operator.gt,
       '>=' : operator.ge}


def _clean_sift_polyphen(vdf, key, to_numeric=False):
    """
    Clean VEP values for entries like SIFT and PolyPhen that provide compound
    classifications and numeric scores as string
    """

    vals = vdf.loc[:, key]

    if to_numeric:
        def __numeric_clean(x):
            if x == '':
                return x
            else:
                return float(re.split('\(|\)', x[1]))
        vdf[key] = vals.apply(__numeric_clean)
    else:
        vdf[key] = vals.apply(lambda x: x.split('(')[0])

    return vdf


def eval_criteria(record, key, criterion, vep_map):
    """
    Checks whether a record meets the criteria specified by key & criterion
    Returns boolean True/False and a list of gene symbols meeting criteria
    """

    # Check for nested keys
    keys = key.split(':')

    # Load relevant operator
    op = ops[criterion[1]]

    # Handle all VEP-related criteria
    if keys[0] == 'CSQ':
        # UTRAnnotator needs to be treated separately
        if keys[1].startswith('UTRAnnotator'):
            criteria_met = False
            genes = []
            for idx, entry in vep2df(record, vep_map, return_dict=True).items():
                if entry[keys[1]] == '':
                    continue
                utra_dict = {v.split(':')[0] : v.split(':')[1] for v in entry[keys[1]].split('&')}
                val = utra_dict.get(keys[2])
                if val is not None:
                    if op(val, criterion[0]):
                        criteria_met = True
                        genes.append(entry['SYMBOL'])

        # Otherwise, the process is the same for all other VEP criteria
        else:
            vdf = vep2df(record, vep_map)
            if keys[1] in 'SIFT PolyPhen'.split():
                vdf = _clean_sift_polyphen(vdf, keys[1], str(criterion[0]).isnumeric())
            query_vals = vdf[keys[1]]
            query_vals = query_vals[(query_vals != '') & (~query_vals.isna())]
            if criterion[1] in '> >= < <='.split():
                query_vals = query_vals.astype(float, errors='ignore')
            hits = query_vals.apply(lambda x: op(x, criterion[0]))
            if hits.any():
                criteria_met = True
                genes = set(vdf.loc[hits.index, 'SYMBOL'].values.tolist())
            else:
                criteria_met = False
                genes = []

    # Handle gene-keyed annotations
    elif key in 'SpliceAI promoter':
        val = record.info.get(key)
        if val is None:
            criteria_met = False
            genes = []
        else:
            criteria_met = True
            genes = [re.split(':|\|', s) for s in val]
            genes = [g for s in genes for g in s]

    # Handle tissue-specific annotations (enhancers and eQTLs)
    elif keys[0] in 'enhancer GTEx_eQTL'.split():
        val = record.info.get(keys[0])
        if val is None:
            criteria_met = False
            genes = []
        else:
            vdf = pd.DataFrame([v.split('|') for v in val])
            if keys[0] == 'enhancer':
                vdf.columns = 'gene score tissue'.split()
                vdf['score'] = vdf.score.astype(float)
            elif keys[0] == 'GTEx_eQTL':
                vdf.columns = 'gene tissue beta'.split()
                vdf['beta'] = vdf.beta.astype(float)
            try:
                hits = vdf[keys[1]].apply(lambda x: op(x, criterion[0]))
            except:
                import pdb; pdb.set_trace()
            if hits.any():
                criteria_met = True
                genes = vdf.gene[hits].values.tolist()
            else:
                criteria_met = False
                genes = []

    # Default behavior: assume single key refers to INFO field
    else:
        val = record.info.get(key)
        genes = []
        if val is None:
            criteria_met = False
        else:
            criteria_met = op(val, criterion[0])

    # Infer which RAS gene for distance criteria based on chromosome
    if key == 'RAS_distance' and criteria_met:
        genes = [{'1' : 'NRAS', 'chr1' : 'NRAS',
                  '11' : 'HRAS', 'chr11' : 'HRAS',
                  '12' : 'KRAS', 'chr12' : 'KRAS'}[record.chrom]]

    return criteria_met, list(set(genes))


def check_single_rule(record, rule_dict, vep_map):
    """
    Check if record meets ALL criteria specified in rule_dict
    Returns boolean True/False and a list of gene symbols meeting criteria
    """

    criteria_res = [eval_criteria(record, key, crit, vep_map) \
                    for key, crit in rule_dict.items()]

    criteria_met = [b for b, g in criteria_res]
    if all(criteria_met):
        glists = [g for b, g in criteria_res]
        genes = list(set([g for l in glists for g in l]))
    else:
        genes = []

    return all(criteria_met), genes


def check_rules(record, rules, vep_map):
    """
    Determine whether a record meets inclusion criteria for a set
    Input rules: must be a list of dicts following the .json spec (see README).
                 Each dict within the list will be evaluated separately
    Returns boolean True/False if any of the dicts in rules are satisfied
    Also returns list of genes meeting any of the rules
    """

    rule_res = [check_single_rule(record, r, vep_map) for r in rules]

    rules_met = [b for b, g in rule_res]
    if any(rules_met):
        glists = [g for b, g in rule_res]
        genes = list(set([g for l in glists for g in l]))
    else:
        genes = []

    return any(rules_met), genes


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='input .vcf')
    parser.add_argument('--sets-json', required=True, help='.json of rules for ' +
                        'each set. See README for specs.')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Load set rules as dict
    with open(args.sets_json) as jfile:
        rules_dict = json.load(jfile)

    # Open connection to input vcf
    if args.vcf in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile('-', 'r')
    else:
        invcf = pysam.VariantFile(args.vcf)

    # Load VEP parser info
    vep_map = parse_vep_map(invcf)

    # Iterate over each record in --vcf and assign to any sets where it meets the criteria
    set_vids = {sid : {} for sid in rules_dict.keys()}
    for record in invcf.fetch():
        for sid, rules in rules_dict.items():
            meets_rules, genes = check_rules(record, rules, vep_map)
            if meets_rules:
                for gene in genes:
                    if gene not in set_vids[sid].keys():
                        set_vids[sid][gene] = set()
                    set_vids[sid][gene].add(record.id)

    # Write results to --outfile
    if args.outfile in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    outfile.write('set_id\tcriteria\tgene\tvids\n')
    for sid, vals in set_vids.items():
        for gene, vids in vals.items():
            if len(vids) > 0:
                vid_str = ','.join(sorted(list(vids)))
                outfile.write('\t'.join([gene + '_' + sid, sid, gene, vid_str]) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()

