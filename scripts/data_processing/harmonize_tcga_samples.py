#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Unify TCGA sample/donors lists between exomes, genotyped arrays, and imputed arrays
"""


import argparse
import pandas as pd
from sys import stdout


def load_exome_samples(exome_ids_in, exome_id_map_in):
    """
    Load list of exome samples and split to donor-indexed dataframe
    """

    import pdb; pdb.set_trace()




def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exome-ids', help='List of exome sample', required=True)
    parser.add_argument('--exome-id-map', help='.tsv mapping exome CRAM/BAM IDs ' +
                        'to their TCGA sample IDs', required=True)
    parser.add_argument('--array-typed-ids', help='List of array-genotyped ' + 
                        'samples', required=True)
    parser.add_argument('--array-imputed-ids', help='List of array-imputed ' + 
                        'samples', required=True)
    parser.add_argument('--tcga-tss-table', help='TCGA TSS code table', required=True)
    parser.add_argument('--out-prefix', help='Prefix for all output files', 
                        required=True)
    args = parser.parse_args()

    # Load exome samples
    ex_df = load_exome_samples(args.exome_ids, args.exome_id_map)


if __name__ == '__main__':
    main()

