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

import logger

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
        if df.empty:
            continue
        df = df[['stringId_A', 'stringId_B', 'preferredName_A', 'preferredName_B', 'score']]
        df = df[df.score >= score_cutoff]
        interactions.append(df)
    if len(interactions) == 0:
        return
    interactions = pd.concat(interactions)
    interactions = interactions.rename(columns={
        'stringId_A': f'id_{level}',
        'stringId_B': f'id_{level + 1}',
        'preferredName_A': f'level{level}',
        'preferredName_B': f'level{level + 1}',
        'score': f'score_{level}'})
    return interactions

def query_string(gene_df, levels, score, output_fp, logger):
    string_version = 'https://string-db.org/api/json/version'
    # TODO: Log STRING version
    ids_df = get_string_id(gene_df['gene'].tolist())
    gene_list = [gene_df['gene'].drop_duplicates().tolist()]
    graph = []
    df = get_string_interaction_partners(ids_df['stringId'].tolist(), score, 0)
    if df is None or df.empty:
        #sys.exit('EXIT No PPIN found for gene list.')
        return gene_list, pd.DataFrame()
    df['token'] = df.apply(lambda row: ''.join(sorted((row['level0'], row['level1']))), axis=1)
    df = df.drop_duplicates('token').drop(columns=['token'])
    gene_list.append(df[~df['level1'].isin(gene_list[0])]['level1'].drop_duplicates().tolist())
    graph.append((df[['level0', 'level1', 'score_0']].drop_duplicates()
                  .rename(columns={'level0': 'geneA',
                                   'level1': 'geneB',
                                   'score_0': 'score'})
                  .assign(level = '0')))
    #logger.write('Querying STRINGdb for interactions between...')
    for level in range(1, levels):
        level_df = get_string_interaction_partners(
            df[f'id_{level}'].drop_duplicates().tolist(), score, level)
        if level_df is None or level_df.empty:
            break
        level_df['token'] = level_df.apply(
            lambda row: ''.join(sorted((row[f'level{level}'], row[f'level{level+1}']))), axis=1)
        level_df = level_df.drop_duplicates('token').drop(columns='token')
        for i in range(level + 1):
            level_df = level_df[~(level_df[f'level{level + 1}'].isin(gene_list[i]))]# == df[f'id_{i}'])]
        graph.append((level_df[[f'level{level}', f'level{level + 1}', f'score_{level}']]
                      .rename(columns={f'level{level}': 'geneA',
                                       f'level{level+1}': 'geneB',
                                       f'score_{level}': 'score'})
                      .assign(level = f'{level}')))
        gene_list.append(level_df[f'level{level + 1}'].drop_duplicates().tolist())
        df = df.merge(level_df, how='right', on=[f'level{level}', f'id_{level}'])
    df = df.drop_duplicates()
    graph = pd.concat(graph)
    for level in range(len(gene_list)):
        write_results(pd.DataFrame({'gene': gene_list[level]}),
                      os.path.join(os.path.dirname(output_fp), f'level{level}_genes.txt'))
    write_results(graph.drop_duplicates(),
                  os.path.join(os.path.dirname(output_fp), 'graph.txt'))
    #write_results(df.drop_duplicates(), output_fp)
    return gene_list, graph


def write_results(df, output_fp):
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    #logger.write('Writing PPIN...')
    df.to_csv(output_fp, sep='\t', index=False)



def parse_input(inputs, logger):
    '''Return a dataframe of gene input.'''
    #logger.write('Parsing input...')
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
    os.makedirs(args.output, exist_ok=True)
    logger = logger.Logger(logfile=os.path.join(args.output, 'query_string.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    gene_df = parse_input(args.input, logger)
    genes, graph = query_string(gene_df, args.levels, args.score, args.output, logger)
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
    
