#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Build table of unique coding consequences present in a VEP-annotated VCF
"""


import argparse
import numpy as np
import pandas as pd
import pysam
from sys import stdin, stdout
from vep_utils import parse_vep_map, vep2df


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input VEP-annotated .vcf')
    parser.add_argument('outfile', help='output .tsv [default: stdout]', default='stdout')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile('-', 'r')
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Prepare metadata for VEP parsing
    vep_map = parse_vep_map(invcf)
    csq_df = pd.DataFrame(columns='set_id gene transcript codon change csq vids'.split())

    # Iterate over records
    for record in invcf.fetch():
        if 'CSQ' not in record.info.keys():
            continue
        
        vdf = vep2df(record, vep_map)
        
        for idx, vals in vdf.iterrows():
            if vals.Feature_type != 'Transcript':
                continue

            gene = vals.SYMBOL
            tx = vals.Feature
            codon = vals.Protein_position
            if codon == '':
                continue
            if vals.HGVSp != '':
                VSp = vals.HGVSp.split(':')[-1]
            else:
                VSp = vals.Consequence
            consequence = vals.Consequence
            set_id = '{}_{}'.format(tx, VSp)

            prev_idx = (csq_df.gene == gene) & \
                       (csq_df.transcript == tx) & \
                       (csq_df.codon == codon) & \
                       (csq_df.change == VSp) & \
                       (csq_df.csq == consequence)
            if prev_idx.any():
                csq_df[prev_idx].vids.map(lambda x: x.add(record.id))
            else:
                csq_df = csq_df.append({'set_id' : set_id, 'gene' : gene, 
                                        'transcript' : tx, 'codon' : codon, 
                                        'change' : VSp, 'csq' : consequence,
                                        'vids' : set([record.id])}, 
                                        ignore_index=True)
    
    # Collapse vid column and write to outfile
    csq_df['vids'] = csq_df.vids.str.join(',')
    try:
        csq_df['codon'] = \
            csq_df.codon.str.split('-').map(lambda x: np.nanmin(np.array([v for v in x if v.isnumeric()]).astype(int)))
    except:
        import pdb; pdb.set_trace()
    csq_df.sort_values(by='gene codon transcript csq'.split()).\
           to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()

