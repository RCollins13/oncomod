#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
General utility functions for RAS modifier project
"""


import pandas as pd


def load_tx_map(tx_tsv):
    """
    Load --transcript-info as a dict
    """
    
    # Load data
    tx_df = pd.read_csv(tx_tsv, sep='\t', header=None)
    tx_df.columns = 'ENST ENSG symbol tx_len'.split()

    # Remove Ensembl version info from ENSG and ENST IDs
    for ecol in 'ENST ENSG'.split():
        tx_df[ecol] = tx_df[ecol].str.replace('.[0-9]+$', '', regex=True)

    # Map everything vs. ENST IDs
    tx_df.set_index('ENST', drop=True, inplace=True)
    
    return {'ENSG' : tx_df.ENSG.to_dict(),
            'ENSG_to_symbol' : tx_df.symbol.set_axis(tx_df.ENSG).to_dict(),
            'symbol' : tx_df.symbol.to_dict(),
            'tx_len' : tx_df.tx_len.to_dict()}

