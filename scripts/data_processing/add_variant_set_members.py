#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Add variant set members to a list of variant sets
"""


import argparse
import pandas as pd
from sys import stdin, stdout


def load_members(memberships_in):
    """
    Load a dict mapping variant set IDs to their constituent variant IDs
    """

    members = {}
    
    for fin in memberships_in:
        map_df = pd.read_csv(fin, sep='\t').iloc[:, [0, -1]]
        map_df.iloc[:, 1] = map_df.iloc[:, 1].apply(lambda x: set(x.split(',')))
        for sid, vids in map_df.itertuples(index=False, name=None):
            if sid in members.keys():
                members[sid].update(vids)
            else:
                members[sid] = vids

    return members


def get_members(setid, members, return_as_string=True):
    """
    Map the constituent variant members onto a single set ID
    """

    if '|' in setid:

        mems = [get_members(subid, members, False) \
                for subid in setid.split('|')]
        mems = [m for m in mems if m is not None]

        if return_as_string:
            return '|'.join([','.join(sorted(list(m))) for m in mems])
    
    else:

        if setid == 'COMUT':
            mems = None

        elif setid in members.keys():
            mems = members[setid]

        else:
            mems = set(list(setid.split(';')))

    if return_as_string:
        return ','.join(sorted(list(mems)))
    else:
        return mems


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set-list', required=True, help='List of variant set ' +
                        'IDs to populate. Will be updated in-place unless stdin ' +
                        'is specified, in which case it will write to stdout.')
    parser.add_argument('--memberships', action='append', help='.tsv mapping ' +
                        'set IDs (first column) to constituent variant IDs ' +
                        '(final column). Can be specified multiple times.')
    args = parser.parse_args()

    # Step 1: read set list into memory
    if args.set_list in 'stdin /dev/stdin -'.split():
        fin = stdin
        fout = stdout
    else:
        fin = open(args.set_list)
        fout = args.set_list
    slist = [l.rstrip() for l in fin.readlines()]

    # Step 2: load memberships into memory
    members = load_members(args.memberships)

    # Step 3: map members onto set IDs
    mappings = {s : get_members(s, members) for s in slist}

    # Step 4: overwrite args.setlist with mappings
    pd.DataFrame.from_dict(mappings, orient='index').\
                 to_csv(fout, sep='\t', index=True, header=False)


if __name__ == '__main__':
    main()

