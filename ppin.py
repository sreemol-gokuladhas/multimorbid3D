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
    if len(ids) == 0:
        return None
    ids = pd.concat(ids)
    try:
        return ids[['queryItem', 'stringId']]
    except:
        return None


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
        try:
            res = requests.post(request_url, data=params)
            time.sleep(10)
            inp = StringIO(res.text.strip(), newline="\n")
            df = pd.read_csv(inp, sep="\t")
            if df.empty:
                continue
            df = df[['stringId_A', 'stringId_B', 'preferredName_A', 'preferredName_B', 'score']]
            df = df[df.score >= score_cutoff]
            interactions.append(df)
        except ConnectionError as e:
            exit(e)
    if len(interactions) == 0:
        return None
    interactions = pd.concat(interactions)
    if interactions.empty:
        return None
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
    return gene_list, graph

def q_string(genes, score, level):
    ids_df = get_string_id(genes)
    if ids_df is None:
        #sys.exit('EXIT No PPIN found for gene list.')
        return pd.DataFrame(columns=['geneA', 'geneB'])
    df = get_string_interaction_partners(ids_df['stringId'].tolist(), score, level)

    if df is None:
        return pd.DataFrame(columns=['geneA', 'geneB'])
    if df.empty:
        return pd.DataFrame(columns=['geneA', 'geneB'])
    try:
        df['token'] = df.apply(
            lambda row: ''.join(sorted((row[f'level{level}'], row[f'level{level+1}']))), axis=1)
    except:
        print(df)
        exit()
    df = (df
          .drop_duplicates('token')
          #.drop(columns=['token', f'id_{level}', f'id_{level+1}'])
          .rename(columns={f'level{level}': 'geneA',
                           f'level{level+1}': 'geneB',
                           f'score_{level}': 'score'})
          .assign(db = 'string')
          )[['geneA', 'geneB', 'score', 'db']]
    return df

def q_proper(genes, proper):
    df1 = proper.loc[proper['Gene1'].isin(genes)]
    #df1 = df1.rename(columns={'Gene1': f'level{level}', 'Gene2': f'level{level + 1}'})
    df1 = df1.rename(columns={'Gene1': f'geneA', 'Gene2': f'geneB'})
    df2 = proper.loc[proper['Gene2'].isin(genes)] 
    df2 = df2.rename(columns={'Gene2': f'geneA', 'Gene1': f'geneB'})
    df = pd.concat([df1, df2]).reset_index(drop=True)
    df = df[df['background_contamination'] == 0]
    if df.empty:
        return pd.DataFrame(columns=['geneA', 'geneB'])
    df['token'] = df.apply(lambda row: ''.join(sorted((row['geneA'], row['geneB']))), axis=1)
    df = (df
          .drop_duplicates('token')
          .drop(columns=['token'])
          .assign(db = 'proper')
          )
    return df

def make_ppin(gene_df, levels, output_dir, ppin_db, score, logger, bootstrap=False):
    res = {}
    graph = []
    gene_list = [gene_df['gene'].drop_duplicates().tolist()]
    for db in ppin_db:
        res[db] = {}
        if db == 'proper':
            proper_fp = os.path.join(os.path.dirname(__file__), 'data/PROPER_v1.csv')
            proper = pd.read_csv(proper_fp).rename(columns={
                'Cell line specificity': f'cell_line',
                'Odds ratio': f'odds_ratio',
                'BH-corrected p-value': f'pval',
                'Read count': 'read_count',
                'Potential background contamination': 'background_contamination'})
            res[db]['db'] = proper
    for level in range(1, levels + 1):
        level_df = []
        for db in res:
            if db == 'string':
                level_df.append(q_string(gene_list[level-1], score, level-1))
            if db == 'proper':
                level_df.append(q_proper(gene_list[level-1], res[db]['db']))
        level_df = pd.concat(level_df)
        if level_df.empty:
            gene_list.append([])
            continue
        level_df['level'] = level
        if not bootstrap:
            graph.append(level_df)
        for i in range(len(gene_list)):
            level_df = level_df[~(level_df[f'geneB'].isin(gene_list[i]))]
        gene_list.append(level_df[f'geneB'].drop_duplicates().tolist())
    if not bootstrap:
        if len(graph) == 0:
            return gene_list, pd.DataFrame(columns=['geneA', 'geneB'])
        graph = pd.concat(graph)
        write_results(graph.drop_duplicates(),
                      os.path.join(output_dir, 'graph.txt'))
        del graph
    for level in range(len(gene_list)):
        write_results(pd.DataFrame({'gene': gene_list[level]}),
                      os.path.join(output_dir, f'level{level}_genes.txt'))
    return gene_list




def query_proper(genes, levels, output_dir, logger):
    proper_fp = os.path.join(os.path.dirname(__file__), 'data/PROPER_v1.csv')
    proper = pd.read_csv(proper_fp).rename(columns={
        'Cell line specificity': f'cell_line',
        'Odds ratio': f'odds_ratio',
        'BH-corrected p-value': f'pval',
        'Read count': 'read_count',
        'Potential background contamination': 'background_contamination'})
    gene_list = [genes['gene'].drop_duplicates().tolist()]
    graph = []
    df_all = pd.DataFrame({'level0': genes['gene'].tolist()})
    #logger.write('Querying PROPER database...')
    for level in range(levels):
        df1 = proper.loc[proper['Gene1'].isin(gene_list[int(level)])]
        #df1 = df1.rename(columns={'Gene1': f'level{level}', 'Gene2': f'level{level + 1}'})
        df1 = df1.rename(columns={'Gene1': f'geneA', 'Gene2': f'geneB'})
        df2 = proper.loc[proper['Gene2'].isin(gene_list[int(level)])] 
        df2 = df2.rename(columns={'Gene2': f'geneA', 'Gene1': f'geneB'})
        df = pd.concat([df1, df2]).reset_index(drop=True)
        df = df[df['background_contamination'] == 0]
        if df.empty:
            break
        df['token'] = df.apply(lambda row: ''.join(sorted((row['geneA'], row['geneB']))), axis=1)
        df = df.drop_duplicates('token').drop(columns=['token'])
        df['level'] = level
        graph.append(df)
        for i in range(level + 1):
            df = df[~(df[f'geneB'].isin(gene_list[i]))]
        gene_list.append(df[f'geneB'].drop_duplicates().tolist())
    if len(graph) == 0:
        return gene_list, pd.DataFrame()
    graph = pd.concat(graph)
    for level in range(len(gene_list)):
        write_results(pd.DataFrame({'gene': gene_list[level]}),
                      os.path.join(os.path.dirname(output_dir), f'level{level}_genes.txt'))
    write_results(graph.drop_duplicates(),
                  os.path.join(os.path.dirname(output_dir), 'graph.txt'))
    
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
        if not 'gene' in df.columns:
            sys.exit('EXITING: No column named "gene" found in input.')
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
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    parser.add_argument(
        '-l', '--levels', default=1, type=int,
        help='Path length (i.e. number of nodes) to query. Default: 1')
    parser.add_argument(
        '-p', '--ppin', nargs='+', default=['string', 'proper'], choices=['string', 'proper'],
        help='PPIN database(s) to query. Default: "all"')
    parser.add_argument(
        '-s', '--score', default=0.7, type=float,
        help='Cut-off score for STRING interactions. Default: 0.7')

    return parser.parse_args()
    
if __name__=='__main__':
    pd.options.mode.chained_assignment = None 
    start_time = time.time()
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    logger = logger.Logger(logfile=os.path.join(args.output, 'ppin.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    gene_df = parse_input(args.input, logger)
    gene_list, graph = make_ppin(gene_df, args.levels, args.output, args.ppin, args.score, logger)

    
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
