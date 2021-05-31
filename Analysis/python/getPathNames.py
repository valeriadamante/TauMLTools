#!/usr/bin/env python
# Get dataId hash
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import ROOT
import pandas
import re

def LoadIdFrames(file_name, id_collections):
    """Load pandas DataFrame from root file."""
    file = ROOT.TFile(file_name, "READ")
    aux = file.Get('L2Summary')
    id_values = {}
    name_values = {}
    id_vector = ROOT.vector('int')()
    for col in id_collections:
        id_values[col] = []
        name_values[col] = []
    for entry in aux:
        for col in id_collections:
            #id_vector += entry #getattr(entry, col + '_hashes')
            name_vector = getattr(entry, col + '_names')
            for n in range(name_vector.size()):
                name = str(name_vector.at(n))
                #if name in name_values[col]:
                #    continue 
                id_values[col].append(n)
                name_values[col].append(name)
    data_frames = {}
    for col in id_collections:
        data_frames[col] = pandas.DataFrame(data = {
            'id' : id_values[col], 'name': name_values[col]
        })
    return data_frames

def GetDataIdHash(df, dataId, use_regex=False):
    """Retrieve hash for the given dataId."""
    id_df = df[df.name.str.contains(dataId, regex=use_regex)]
    result = {}
    for n in range(id_df.shape[0]):
        result[id_df.name.values[n]] = id_df.id.values[n]
    return result

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get dataId hash.')
    parser.add_argument('input', type=str, help="Input file with anaTuple")
    parser.add_argument('dataId', type=str, help="dataId")
    args = parser.parse_args()

    id_collections = ['module']
    data_frames = LoadIdFrames(args.input, id_collections)
    id_found = False
    for col in id_collections:
        df = data_frames[col]
        ids = GetDataIdHash(df, args.dataId, False)
        if len(ids) == 0:
            ids = GetDataIdHash(df, args.dataId, True)
        if len(ids) != 0:
            id_found = True
            print('{} IDs:'.format(col))
            for name in ids:
                print('\t{:<20}\t{}'.format(ids[name], name))
    if not id_found:
        print("ERROR: Id not found")
