#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Reformat a list of mutations as a .json to be used by define_control_samples.py
"""


import argparse
import json
import re
from string import ascii_uppercase


def reformat_mut(mut):
    """
    Reformat a single mutation consequence string as a dict of criteria to apply
    """

    # Define protein position
    aa_pos = int(re.split('[A-z_]', mut)[1])
    mdict = {'Protein_position' : aa_pos}

    # Other annotations depend on variant type
    alt = re.split('[0-9]+', mut)[1]
    if 'del' == alt.lower():
        mdict['variant_type'] = 'deletion'
    elif 'ins' == alt.lower():
        mdict['variant_type'] = 'insertion'
    elif 'delins' == alt.lower():
        mdict['variant_type'] = 'delins'
    elif alt == '*':
        mdict['Consequence'] = 'stop_gained'
    elif alt in ascii_uppercase:
        mdict['Consequence'] = 'missense_variant'
        mdict['Amino_acids'] = alt

    return mdict


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mutations', help='input .txt file of mutations')
    parser.add_argument('outfile', help='Output .json')
    parser.add_argument('--gene', default='KRAS', help='Gene corresponding to ' +
                        '--mutations')
    args = parser.parse_args()

    out_dict = {args.gene : []}

    # Iterate over mutations
    with open(args.mutations) as fin:
        for mut in fin.readlines():
            mut = re.sub('^p\.', '', mut.rstrip())

            # Reformat each mutation as .json-style criteria
            mut_fmt = reformat_mut(mut)
            out_dict[args.gene].append(mut_fmt)

    # Dump as .json to outfile
    with open(args.outfile, 'w') as fout:
        json.dump(out_dict, fout, ensure_ascii=False)


if __name__ == '__main__':
    main()

