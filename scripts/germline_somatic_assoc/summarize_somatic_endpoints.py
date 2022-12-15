#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Summarize number of somatic endpoints to test for each cancer type & gene
"""


import argparse
import pandas as pd
from sys import stdout


cancers = 'PDAC CRAD SKCM LUAD'.split()
ras_genes = 'NRAS HRAS KRAS'.split()
category_descriptions = \
    {'mutations' : '1. Frequent RAS mutations',
     'codons'    : '2. Recurrently mutated RAS codons',
     'burden'    : '3. Collapsed mutation sets'}


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutations', action='append', nargs='+', 
                        help='frequencies corresponding to individual mutations')
    parser.add_argument('--codons', action='append', nargs='+', 
                        help='frequencies corresponding to recurrent codon mutations')
    parser.add_argument('--burden-sets', action='append', nargs='+', 
                        help='frequencies corresponding to burden sets')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Build dict for collecting results
    res = {cat : {cncr : {gene : set() for gene in ras_genes} for cncr in cancers} \
           for cat in category_descriptions.keys()}
    
    # Load individual mutations
    for infile in args.mutations:
        res['mutations'] = update_res(res['mutations'], infile)





if __name__ == '__main__':
    main()

