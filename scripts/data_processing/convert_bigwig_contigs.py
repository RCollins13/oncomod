#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Convert a bigWig from UCSC-style ("chr") to Ensembl-style (no "chr") contigs
"""


import pandas as pd
import pyBigWig as pbw
import sys


# Open connections to in/out files
infile, outfile = sys.argv[1:3]
old = pbw.open(infile)
new = pbw.open(outfile, 'w')

# Add header to outfile
contig_map = {k : k.replace('chr', '') for k in old.chroms().keys()}
new.addHeader([(contig_map[k], v) for k, v in old.chroms().items()])

# Process each chromosome in chunks
step_bp = 1000000
for contig, n_bp in old.chroms().items():
    for start in range(0, n_bp, step_bp):
        end = min([start + step_bp, n_bp])
        vals = pd.DataFrame(old.intervals(contig, start, end),
                            columns='start end value'.split())
        vals['contig'] = contig_map[contig]
        new.addEntries(vals.contig.values.tolist(),
                       vals.start.values.tolist(),
                       ends=vals.end.values.tolist(),
                       values=vals.value.values.tolist())

# Close new to clear buffer
new.close()