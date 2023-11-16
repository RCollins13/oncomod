#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Clean up verbose VEP output for OncoMod VCFs
"""


import argparse
import pybedtools as pbt
from sys import stdout


def dump_cache(tx_dict, outfile):
    """
    Write all currently stored info to outfile
    """

    # Compute transcript length for each transcript and write to outfile
    for txid, txvals in tx_dict.items():
        if txvals.get('ENSG') is None:
            continue
        ex_bt = pbt.BedTool('\n'.join(tx_dict[txid]['tx_str']), from_string=True)
        tx_len = sum(map(len, ex_bt.sort().merge()))
        outvals = [str(x) for x in [txid, txvals['ENSG'], txvals['symbol'], tx_len]]
        outfile.write('\t'.join(outvals) + '\n')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', help='input gtf', required=True)
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    parser.add_argument('--no-header', action='store_true', help='no header on output')
    args = parser.parse_args()

    # Open connection to outfile
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')
    if not args.no_header:
        outfile.write('#ENST\tENSG\tsymbol\ttx_length\n')

    # Load info from transcripts and exons into memory
    prev_contig = ''
    tx_dict = {}
    for feature in pbt.BedTool(args.gtf).\
                       filter(lambda x: x.fields[2] in 'transcript exon'.split()):

        if feature.chrom != prev_contig:
            dump_cache(tx_dict, outfile)
            tx_dict = {}
        prev_contig = feature.chrom
        
        ftype = feature.fields[2]
        txid = feature.attrs.get('transcript_id')
        if txid is None:
            continue

        if txid not in tx_dict.keys():
            tx_dict[txid] = {'tx_str' : []}
        
        if ftype == 'transcript':
            tx_dict[txid].update({'ENSG' : feature.attrs.get('gene_id'),
                                  'symbol' : feature.attrs.get('gene_name')})

        if ftype == 'exon':
            cvals = [feature.chrom, feature.start, feature.end]
            tx_dict[txid]['tx_str'].append('\t'.join([str(x) for x in cvals]))

    # Close outfile to clear buffer
    dump_cache(tx_dict, outfile)
    outfile.close()


if __name__ == '__main__':
    main()

