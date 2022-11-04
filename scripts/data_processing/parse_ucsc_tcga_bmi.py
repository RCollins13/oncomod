#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Parse TCGA BMI data streamed from UCSC
Assumes stdin and stdout for I/O
"""

from sys import stdin, stdout

for line in stdin.readlines():
    samples, values = line.rstrip().split('\t')
    for sample, value in zip(samples.split(','), values.split(',')):
        if value != '--':
            stdout.write('{}\t{}\n'.format(sample, value))

