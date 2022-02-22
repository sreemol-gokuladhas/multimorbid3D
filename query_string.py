#! /usr/bin/env python

import pandas as pd
#import stringdb
import argparse
import requests
import sys
import os
import stringdb_params as sdb
from io import StringIO
import time


def get_string_id(genes):
    ''' Get STRING IDs of a list of gene
    Returns a pandas DataFrame with columns, queryItem and stringId
    '''
    ids = []
    chunksize = 2000
    gene_chunks = [genes[x:x+chunksize] for x in range(0, len(genes), chunksize)]
    for chunk in gene_chunks:
        params = sdb.get_params(genes)
        request_url = sdb.get_request_ids_url()
        res = requests.post(request_url, data=params)
        time.sleep(1)
        inp = StringIO(res.text.strip(), newline="\n")
        df = pd.read_csv(inp, sep="\t")
        ids.append(df)
    ids = pd.concat(ids)
    return ids[['queryItem', 'stringId']]


def get_string_interaction_partners(genes, score_cutoff, level):
    ''' Get STRING IDs of a list of gene
    Returns a pandas DataFrame with columns, queryItem and stringId
    '''
    interactions = []
    chunksize = 1000
    gene_chunks = [genes[x:x+chunksize] for x in range(0, len(genes), chunksize)]
    for chunk in gene_chunks:
        params = sdb.get_params(chunk)
        request_url = sdb.get_request_interaction_url()
        res = requests.post(request_url, data=params)
        time.sleep(1)
        inp = StringIO(res.text.strip(), newline="\n")
        df = pd.read_csv(inp, sep="\t")
        df = df[['stringId_A', 'stringId_B', 'preferredName_A', 'preferredName_B', 'score']]
        df = df[df.score >= score_cutoff]
        interactions.append(df)
    interactions = pd.concat(interactions)
    interactions = interactions.rename(columns={
        'stringId_A': f'id_{level}',
        'stringId_B': f'id_{level + 1}',
        'preferredName_A': f'level{level}',
        'preferredName_B': f'level{level + 1}',
        'score': f'score_{level}'})
    return interactions

def query_string(gene_df, levels, score, output_fp):
    ids_df = get_string_id(gene_df['gene'].tolist())
    df = get_string_interaction_partners(ids_df['stringId'].tolist(), score, 0)
    #print('Querying STRINGdb for interactions between...')
    for level in range(1, levels):
        #print(f'\tLevels {level-1} and {level}')
        level_df = get_string_interaction_partners(
            df[f'id_{level}'].drop_duplicates().tolist(), score, level)
        df = df.merge(level_df, how='right', on=[f'level{level}', f'id_{level}'])
        for i in range(level):
            df = df[~(df[f'id_{level + 1}'] == df[f'id_{i}'])]
    df = df.drop_duplicates()
    write_results(df, output_fp)
    return df


def write_results(df, output_fp):
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    #print('Writing PPIN...')
    df.to_csv(output_fp, sep='\t', index=False)



def parse_input(inputs):
    '''Return a dataframe of gene input.'''
    #print('Parsing input...')
    df = pd.DataFrame()
    if os.path.isfile(inputs[0]): # Input is file.
        df = pd.read_csv(inputs[0], sep='\t')
        df = df[['gene']].drop_duplicates()
    else: # Input probably space-separated list of genes.
        df = pd.DataFrame({'gene': [i.upper() for i in inputs]})
        df = df[['gene']].drop_duplicates()
    return df


def parse_args():
    parser = argparse.ArgumentParser(
        description='Query STRINGdb for proteins that interact with input proteins.')
    parser.add_argument(
        '-i', '--input', required=True, nargs='+',
        help='''A space-separated list of gene symbols or filepath to a file 
        containing gene symbols in the 'gene' column.''' )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Filepath to write results.')
    parser.add_argument(
        '-l', '--levels', default=1, type=int,
        help='Path length (i.e. number of nodes) to query. Default: 1')
    parser.add_argument(
        '-s', '--score', default=0.7, type=float,
        help='Cut-off score for STRING interactions. Default: 0.7')

    return parser.parse_args()

    
if __name__=='__main__':
    pd.options.mode.chained_assignment = None 
    start_time = time.time()
    args = parse_args()
    gene_df = parse_input(args.input)
    df = query_string(gene_df, args.levels, args.score, args.output)

    #print('Done.')
    #print(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
