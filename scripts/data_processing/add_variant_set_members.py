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


def get_members(setid, members, return_as_string=True):
    """
    Map the constituent variant members onto a single set ID
    """

    mset = set([])

    for subid in setid.split('|'):
        
        if subid == 'COMUT':
            continue

        elif subid in members.keys():
            mset.update(members[subid])

        else:
            mset.update(set(list(subid.split(';'))))

    if return_as_string:
        return ','.join(sorted(list(mset)))
    else:
        return mset


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set-list', required=True, help='List of variant set ' +
                        'IDs to populate. Will be updated in-place.')
    parser.add_argument('--memberships', action='append', help='.tsv mapping ' +
                        'set IDs (first column) to constituent variant IDs ' +
                        '(final column). Can be specified multiple times.')
    args = parser.parse_args()

    # Step 1: read set list into memory
    with open(args.set_list) as fin:
        slist = [l.rstrip() for l in fin.readlines()]

    # Step 2: load memberships into memory
    members = {}
    for fin in args.memberships:
        map_df = pd.read_csv(fin, sep='\t').iloc[:, [0, -1]]
        map_df.iloc[:, 1] = map_df.iloc[:, 1].apply(lambda x: set(x.split(',')))
        for sid, vids in map_df.itertuples(index=False, name=None):
            if sid in members.keys():
                members[sid].update(vids)
            else:
                members[sid] = vids

    # Step 3: map members onto set IDs
    mappings = {s : get_members(s, members) for s in slist}

    # Step 4: overwrite args.setlist with mappings
    pd.DataFrame.from_dict(mappings, orient='index').\
                 to_csv(args.set_list, sep='\t', index=True, header=False)


if __name__ == '__main__':
    main()

