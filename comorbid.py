#! /usr/bin/env python

import pandas as pd
import istarmap
import argparse
import requests
import sys
import os
import stringdb_params as sdb
from io import StringIO
import time
import shutil
import multiprocessing as mp
from tqdm import tqdm
from itertools import repeat

import query_grn
import query_string
import query_proper
import get_ppi_eqtls
import find_snp_disease



def write_results(df, output_fp):
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    print('Writing output...')
    df.to_csv(output_fp, sep='\t', index=False)



def parse_input(inputs):
    '''Return a dataframe of gene input.'''
    print('Parsing input...')
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
        '-g', '--genes', nargs='+',
        help='''A space-separated list of gene symbols or filepath to a file 
        containing gene symbols in the 'gene' column.''' )
    parser.add_argument(
        '-s', '--snps', nargs='+',
        help='''A space-separated list of SNP rsIDs or filepath to a file 
        containing SNP rsids in the 'snp' column.''' )
    parser.add_argument(
        '--trait',
        help='''GWAS trait to query. Note: this flag is mutually exclusive with
        the --snps and --pmid flags''')
    parser.add_argument(
        '--pmid',
        help='''PubMed ID of the GWAS to query. Note: this flag is mutually exclusive with 
        the --snps and --trait flag''')
    parser.add_argument(
        '--grn-dir', required=True,
        help='''Directory containing tissue gene regulatory network.
        The subdirectories should contain significant_eqtls.txt for each chromosome.''')
    parser.add_argument(
        '--gwas', default=None,
        help='''Filepath to GWAS associations.
        Default: Associations from the GWAS Catalog
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    parser.add_argument(
        '--non-spatial', action='store_true', default=False,
        help='Include non-spatial eQTLs.')
    parser.add_argument(
        '--non-spatial-dir', default='data/GTEx/', help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '-l', '--levels', default=1, type=int,
        help='Path length (i.e. number of nodes) to query. Default: 1')
    parser.add_argument(
        '--ppin', required=True, choices=['string', 'proper'],
        help='''The protein-protein-interaction data to use.''')
    parser.add_argument(
        '--string-score', default=0.7, type=float,
        help='Cut-off score for STRING interactions. Default: 0.7')
    parser.add_argument(
        '--bootstraps', default=1000, type=int,
        help='Number of bootstrap datasets. Default: 1000')

    return parser.parse_args()

def parse_snps(snp_arg, trait_arg, pmid_arg, gwas, grn, output_dir):
    if (snp_arg and trait_arg) or (snp_arg and pmid_arg) or (trait_arg and pmid_arg):
        sys.exit('Only one of --snps, --trait, or --pmid is required.\nExiting.')
    snps = pd.DataFrame()
    print('Parsing SNP input...')
    if snp_arg:
        if os.path.isfile(snp_arg[0]):
            df = pd.read_csv(snp_arg[0], sep='\t')
            snps = df['snp'].drop_duplicates()
        else:
            df = pd.DataFrame({'snp':  snp_arg})
            snps = df['snp'].drop_duplicates()
    elif trait_arg:
        snps = query_grn.extract_trait_snps(trait_arg, gwas)
    elif pmid_arg:
        snps = query_grn.extract_pmid_snps(pmid_arg, gwas)
    eqtls = query_grn.get_eqtls(snps, grn, output_dir)
    
    return snps, eqtls

def parse_genes(genes_args):
    print('Parsing gene input...')
    df = pd.DataFrame()
    if os.path.isfile(genes_args[0]):
        df = pd.read_csv(genes_args[0], sep='\t')
        df = df[['gene']].drop_duplicates()
    else:
        df = pd.DataFrame({'gene': [i.upper() for i in genes_args]})
        df = df[['gene']].drop_duplicates()
    return df

def join_path(*args):
    fp = ''
    for arg in args:
        fp = os.path.join(fp, arg)
    return fp

def pipeline(genes, gwas, output_dir, args):
    # PPIN
    ppin = []
    if args.ppin == 'string':
        ppin = query_string.query_string(
            genes, args.levels, args.string_score, join_path(output_dir, 'ppin.txt'))
    else:
        ppin = query_proper.query_proper(
            genes, args.levels, join_path(output_dir, 'ppin.txt'))

    # PPIN eQTLs
    ppin_eqtls = get_ppi_eqtls.get_snps(ppin, args.grn_dir, output_dir)
    
    # Traits
    sig_res = find_snp_disease.find_disease(gwas, output_dir, output_dir)
    return sig_res

def prep_bootstrap(sim, gene_num, sims_dir, res_dict, grn_genes, gwas, args):
    sim_output_dir = join_path(sims_dir, sim)
    sim_genes = pd.DataFrame(
        {'gene': grn_genes.sample(gene_num, random_state=int(sim)).tolist()})
    sim_res = pipeline(sim_genes, gwas, sim_output_dir, args)
    if sim_res is None:
        return
    for i, row in sim_res.iterrows():
        k = row['level'] + '__' + row['trait']
        try:
            res_dict[k] += 1
        except:
            pass
            
def bootstrap_genes(sig_res, genes, gwas, num_sims, grn, args):
    gene_num = genes[genes['gene'].isin(grn['gene'])]['gene'].nunique()
    sims_dir = join_path(args.output_dir, 'bootstrap')
    manager = mp.Manager()
    res_dict = manager.dict()
    for _, row in sig_res.iterrows():
        k = row['level'] + '__' + row['trait']
        if not k in res_dict.keys():
            res_dict[k] = 0
    desc = 'Boostrapping'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
    sims = [str(i) for i in range(num_sims)]
    with mp.Pool(16) as pool:
        for _  in tqdm(
                pool.istarmap(
                    prep_bootstrap,
                    zip(
                        sims,
                        repeat(gene_num),
                        repeat(sims_dir),
                        repeat(res_dict),
                        repeat(grn['gene'].drop_duplicates()),
                        repeat(gwas),
                        repeat(args)
                    )
                ),
                total=len(sims), desc=desc, bar_format=bar_format,
                unit='simulations', ncols=80
        ):
            pass
    sim_df = []
    for k in res_dict:
        ks = k.split('__')
        sim_df.append([ks[0], ks[1], res_dict[k]])
    sim_df = pd.DataFrame(sim_df, columns=['level', 'trait', 'sim_count'])
    # Using (count + 1) / (num_sims +1) See https://doi.org/10.1086/341527
    sim_df['sim_pval'] = (sim_df['sim_count'] + 1) / (num_sims + 1)
    res = (sig_res
           .merge(sim_df, on=['level', 'trait'], how='left')
           .fillna(0)
    )
    shutil.rmtree(sims_dir)
    write_results(res, join_path(args.output_dir, 'significant_enrichment_bootstrap.txt'))
    
        
if __name__=='__main__':
    pd.options.mode.chained_assignment = None
    args = parse_args()
    if not args.genes and not args.snps and not args.trait and not args.pmid:
        sys.exit('FATAL: One of --genes, --snps, --trait, or --pmid is required.\nExiting.')
    start_time = time.time()

    if args.genes and (args.snps or args.trait or args.pmid):
        sys.exit('Only one of --genes, --snps, --trait, or --pmid is required.\nExiting.')
    gwas = query_grn.parse_gwas(args.gwas)
    grn = query_grn.parse_grn(args.grn_dir)
    snps = []
    genes = []
    if args.genes:
        genes = parse_genes(args.genes)
    else:
        snps, genes = parse_snps(args.snps, args.trait, args.pmid, gwas, grn, args.output_dir)

    sig_res = pipeline(genes, gwas, args.output_dir,  args)
    
    # Bootstrap
    #if args.genes:
    bootstrap_genes(sig_res, genes, gwas, args.bootstraps, grn, args)
    print('Done.')
    print(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
