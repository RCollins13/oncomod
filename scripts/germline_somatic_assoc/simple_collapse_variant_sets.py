#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Simple helper to collapse variant sets from a two-column .tsv input
"""


import argparse
import csv
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tsv-in', default='stdin', metavar='file', 
                        help='two-column .tsv mapping set ID to variant IDs ' + 
                        '[default: stdin]')
    parser.add_argument('--tsv-out', help='output .tsv [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Make dictionary keyed on set ID for collecting results
    res = {}

    # Open connection to input file
    if args.tsv_in in 'stdin /dev/stdin -'.split():
        infile = csv.reader(stdin, delimiter='\t')
    else:
        infile = csv.reader(open(args.tsv_in), delimiter='\t')

    # Iterate over lines in input file
    for sid, vids_str in infile:
        vids = vids_str.split(',')
        if sid not in res.keys():
            res[sid] = set()
        res[sid].update(set(vids))

    # Once finished parsing input file, write to --outfile
    if args.tsv_out in 'stdout /dev/stdout -'.split():
        outfile = stdout
    else:
        outfile = open(args.tsv_out, 'w')
    for sid, vids in res.items():
        outfile.write('{}\t{}\n'.format(sid, ','.join(sorted(list(vids)))))

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

