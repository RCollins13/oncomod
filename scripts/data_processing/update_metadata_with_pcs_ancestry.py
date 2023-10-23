#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Update a cohort's sample metadata manifest with new PCs and ancestry labels
"""


import argparse
import pandas as pd


def update_pops(md, pop_labels_in):
    """
    Update sample ancestry assignments with population labels in pop_labels_in
    """

    # Format labels
    labs = pd.read_csv(pop_labels_in, sep='\t', header=None)
    labs.columns = [md.columns[0], 'POPULATION']
    labs.set_index(labs.columns[0], drop=True, inplace=True)
    labs = labs.POPULATION.to_dict()

    # Map labels onto metadata
    md.loc[:, 'POPULATION'] = md.iloc[:, 0].map(labs)

    return md


def update_pcs(md, pcs_in):
    """
    Update sample PCs based on values in pcs_in
    """

    # Format PCs
    pcs = pd.read_csv(pcs_in, sep='\t')
    pcs.set_index(pcs.columns[0], drop=True, inplace=True)
    pcs = pcs.to_dict()

    # Map PCs onto metadata
    for colname in [k for k in md.columns if k.startswith('PC')]:
        md.loc[:, colname] = md.iloc[:, 0].map(pcs[colname])
    
    return md


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metadata', metavar='.tsv', help='Old sample metadata.tsv',
                        required=True)
    parser.add_argument('--pop-labels', metavar='.tsv', help='Two-column .tsv ' +
                        'mapping sample IDs to population labels', required=True)
    parser.add_argument('--PCs', metavar='.tsv', help='New PCs for all samples', 
                        required=True)
    parser.add_argument('--outfile', help='Path to output file. If not specified, ' +
                        'will overwrite --metadata in place.')
    args = parser.parse_args()

    # Load old metadata
    md = pd.read_csv(args.metadata, sep='\t')

    # Update population labels
    md = update_pops(md, args.pop_labels)

    # Update PC values
    md = update_pcs(md, args.PCs)

    # Write updated metadata to --outfile
    if args.outfile is None:
        outfile = args.metadata
    else:
        outfile = args.outfile
    md.to_csv(args.outfile, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    main()

