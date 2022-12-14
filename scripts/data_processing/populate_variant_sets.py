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
import pandas as pd
import pysam
from sys import stdin, stdout
from vep_utils import parse_vep_map, vep2df


ops = {'==' : operator.eq,
       '!=' : operator.ne,
       '<' : operator.lt,
       '<=' : operator.le,
       '>' : operator.gt,
       '>=' : operator.ge}


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
        vdf = vep2df(record, vep_map)
        hits = vdf[keys[1]].apply(lambda x: op(x, criterion[0]))
        if hits.any():
            criteria_met = True
            genes = set(vdf.SYMBOL[hits].values.tolist())
        else:
            criteria_met = False
            genes = []

    elif len(keys) > 1:
        import pdb; pdb.set_trace()

    # Handle SpliceAI annotations
    elif key == 'SpliceAI':
        val = record.info.get(key)
        if val is None:
            criteria_met = False
            genes = []
        else:
            criteria_met = True
            if len(val) > 1:
                import pdb; pdb.set_trace()
            genes = [s.split(':')[0] for s in val]

    # Default behavior: assume single key refers to INFO field
    else:
        val = record.info.get(key)
        genes = []
        if val is None:
            criteria_met = False
        else:
            criteria_met = op(val, criterion[0])

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
    outfile.write('set_id\tgene\tvids\n')
    for sid, vals in set_vids.items():
        for gene, vids in vals.items():
            if len(vids) > 0:
                vid_str = ','.join(sort(list(vids)))
                outfile.write('\t'.join(sid, gene, vid_str) + '\n')
    outfile.close()


if __name__ == '__main__':
    main()

