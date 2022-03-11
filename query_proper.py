#!/usr/bin/env python
import pandas as pd
import sys
import os
import argparse

def parse_eqtls(eqtls_fp):
    #print('Parsing input...')
    if os.path.isfile(eqtls_fp[0]): 
        df = pd.read_csv(eqtls_fp[0], sep='\t')
        return df['gene'].drop_duplicates()
    else:
        return pd.DataFrame({'gene': eqtls_fp})

def query_proper(genes, levels, output_fp):
    proper_fp = 'data/PROPER_v1.csv'
    proper = pd.read_csv(proper_fp).rename(columns={
            'Cell line specificity': f'cell_line',
            'Odds ratio': f'odds_ratio',
            'BH-corrected p-value': f'pval',
        'Read count': 'read_count',
        'Potential background contamination': 'background_contamination'})
    gene_list = [genes['gene'].drop_duplicates().tolist()]
    graph = []
    df_all = pd.DataFrame({'level0': genes['gene'].tolist()})
    #print('Querying PROPER database...')
    for level in range(levels):
        #print(f'\tlevel {level + 1}...')
        df1 = proper.loc[proper['Gene1'].isin(gene_list[int(level)])]
        #df1 = df1.rename(columns={'Gene1': f'level{level}', 'Gene2': f'level{level + 1}'})
        df1 = df1.rename(columns={'Gene1': f'geneA', 'Gene2': f'geneB'})
        df2 = proper.loc[proper['Gene2'].isin(gene_list[int(level)])] 
        #df2 = df2.rename(columns={'Gene2': f'level{level}', 'Gene1': f'level{level + 1}'})
        df2 = df2.rename(columns={'Gene2': f'geneA', 'Gene1': f'geneB'})
        df = pd.concat([df1, df2]).reset_index(drop=True)
        df = df[df['background_contamination'] == 0]
        if df.empty:
            break
        df['token'] = df.apply(lambda row: ''.join(sorted((row['geneA'], row['geneB']))), axis=1)
        df = df.drop_duplicates('token').drop(columns=['token'])
        df['level'] = level
        '''
        df = df.rename(columns={
            'Cell line specificity': f'cell_line{level}',
            'Odds ratio': f'odds_ratio{level}',
            'BH-corrected p-value': f'pval{level}'})
        '''
        #cols = [f'level{level}', f'level{level + 1}', f'cell_line{level}']#, f'odds_ratio{level}', f'pval{level}']
        #df = df[cols]
        graph.append(df)
        for i in range(level + 1):
            df = df[~(df[f'geneB'].isin(gene_list[i]))]
        gene_list.append(df[f'geneB'].drop_duplicates().tolist())

        #df_all = df_all.merge(df[cols], how='left', on=[f'level{level}'])
    #df_all =  df_all.dropna(subset=['level1']).fillna('')
    for level in range(len(gene_list)):
        write_results(pd.DataFrame({'gene': gene_list[level]}),
                      os.path.join(os.path.dirname(output_fp), f'level{level}_genes.txt'))
    if len(graph) == 0:
        #sys.exit('EXIT: No PPIN found for gene list.')
        return gene_list, pd.DataFrame()
    graph = pd.concat(graph)
    write_results(graph, os.path.join(os.path.dirname(output_fp), 'graph.txt'))
    #write_results(df_all, output_fp)
    return gene_list, graph

def write_results(df, out):
    #print('Writing PPIN...')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep='\t', index=False)
    
def parse_args():
    parser = argparse.ArgumentParser(
        description='Query PROPER-seq for proteins that interact with input proteins.')
    parser.add_argument(
        '-g', '--genes', required=True, nargs='+',
        help='''A space-separated list of gene symbols or filepath to a file 
        containing gene symbols in the 'gene' column.''' )
    parser.add_argument(
        '-o', '--output', required=True, help='Filepath to write results.')
    parser.add_argument('-l', '--levels', default=1, type=int,
                        help='Path length (i.e. number of nodes) to query. Default: 1')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    genes = parse_eqtls(args.input)
    gene_list, graph = query_proper(genes, args.levels, args.output)

