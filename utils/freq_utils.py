#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Utility functions for mutation & mutation set frequency calculations
"""


import pandas as pd


def load_variant_sets(tsv_in):
    """
    Build a dict mapping variant set ID : list of variant IDs in set
    """

    vdf = pd.read_csv(tsv_in, sep='\t')
    vdf.set_index(vdf.columns[0], drop=True, inplace=True)

    return vdf.iloc[:, -1].str.split(',').to_dict()


def load_cancer_map(meta_in):
    """
    Build a dict mapping sample ID to cancer type
    """

    mdf = pd.read_csv(meta_in, sep='\t')
    mdf.set_index(mdf.columns[0], drop=True, inplace=True)

    return mdf['CANCER_TYPE'].to_dict()


def load_dosage_matrix(tsv_in, elig_samples=None):
    """
    Load dosage matrix as pd.DataFrame and subset to elig_samples (if not None)
    """

    # Load matrix and set index as first column
    ad_df = pd.read_csv(tsv_in, sep='\t')
    ad_df = ad_df.set_index(ad_df.columns[0], drop=True).astype(int, errors='ignore')

    # Subset to eligible samples if optioned
    if elig_samples is not None:
        ad_df = ad_df.loc[:, ad_df.columns.isin(elig_samples)]

    return ad_df
