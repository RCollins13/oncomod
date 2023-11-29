#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Compute median value per sample from a two-column tsv of sample, value
"""


import argparse
from numpy import nanmedian
from sys import stdout


def is_gz_file(filepath):
    """
    Check whether input file is gzipped
    Taken from https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    """
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsv_in', help='input .tsv of sample, value observations')
    args = parser.parse_args()

    # Create dictionary for collecting results
    res = {}

    # Open connection to tsv
    if is_gz_file(args.tsv_in):
        import gzip
        fin = gzip.open(args.tsv_in, 'rt')
    else:
        fin = open(args.tsv_in)

    # Iterate over lines in file, adding sample:value pairs to results dictionary
    for line in fin.readlines():
        sample, value = line.rstrip().split('\t')
        if sample not in res.keys():
            res[sample] = []
        res[sample].append(float(value))

    # Once finished, compute median for each sample and write to stdout
    for sample, vals in res.items():
        med = nanmedian(vals)
        stdout.write('{}\t{}\n'.format(sample, med))

    # Clear buffer before exiting
    stdout.close()


if __name__ == '__main__':
    main()

