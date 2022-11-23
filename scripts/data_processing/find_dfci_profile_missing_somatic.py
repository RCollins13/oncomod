#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Identify DFCI-PROFILE donors who are missing somatic mutation data
"""


import argparse
import pandas as pd
from preprocess_dfci_profile_somatic import load_id_map


def load_and_subset_ids(fin, id_map, id_column='SAMPLE_ACCESSION_NBR'):
    """
    Extracts a set of IDs present in infile (fin)
    Converts those IDs to PBP IDs
    """

    vals = pd.read_csv(fin, header=0, low_memory=False)[id_column].map(id_map)

    return set(vals[~vals.isna()].tolist())


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mutation-csv', help='Somatic mutation .csv', required=True)
    parser.add_argument('--cna-csv', help='Copy number alteration .csv', required=True)
    parser.add_argument('--samples-list', help='List of samples to evaluate', required=True)
    parser.add_argument('--id-map-tsv', help='.tsv mapping between various IDs', required=True)
    parser.add_argument('-o', '--out-prefix', help='Prefix for output files', required=True)
    args = parser.parse_args()

    # Load ID map for samples of interest
    id_map = load_id_map(args.samples_list, args.id_map_tsv)
    all_ids = set(id_map.values())

    # Load samples present in both mutation and CNA files
    mut_ids = load_and_subset_ids(args.mutation_csv, id_map)
    cna_ids = load_and_subset_ids(args.cna_csv, id_map)
    somatic_ids = mut_ids.intersection(cna_ids)

    # Find exclusion set of samples lacking both mutation and CNA data
    no_mut = all_ids.difference(mut_ids)
    no_cna = all_ids.difference(cna_ids)
    no_somatic = all_ids.difference(somatic_ids)

    # Write list of samples missing somatic data to output files
    for slist, suffix in [(no_mut, '.mut.list'), (no_cna, '.cna.list'), (no_somatic, '.list')]:
        fout = open(args.out_prefix + suffix, 'w')
        for sample in slist:
            fout.write(sample + '\n')
        fout.close()

if __name__ == '__main__':
    main()

