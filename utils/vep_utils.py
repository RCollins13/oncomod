#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Helper functions for parsing VEP annotations in VCF
"""


import pandas as pd


def parse_vep_map(vcf):
    """
    Parse VEP field mappings to variable names
    Input must be a pysam.VariantFile
    """

    vep_text = vcf.header.info.get('CSQ').description
    vep_fields = vep_text.split('Format: ')[1].replace('\'', '').split('|')

    return {i : k for i, k in enumerate(vep_fields)}


def vep2df(record, vep_map, vep_pop=[], return_dict=False):
    """
    Build pd.DataFrame of all VEP entries in a VCF record (excluding keys in vep_pop)
    """

    vep_vals = {}

    for i, vep_str in enumerate(record.info.get('CSQ', [])):
        vep_vals[i] = {k : v for k, v in zip(vep_map.values(), vep_str.split('|')) \
                       if k not in vep_pop}
    
    if return_dict:
        return vep_vals
    else:
        return pd.DataFrame.from_dict(vep_vals, orient='index')

