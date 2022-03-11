#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm

import ld_proxy


def parse_gwas(gwas_fp):
    #print('Parsing GWAS associations...')
    cols = ['SNPS', 'SNP_ID_CURRENT', 'DISEASE/TRAIT', #'MAPPED_GENE',
            'P-VALUE',  'OR or BETA', '95% CI (TEXT)']
    if gwas_fp is None:
        fp = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'
        df = pd.read_csv(fp, sep='\t', usecols=cols, low_memory=False)
    else:
        if not os.path.isfile(gwas_fp):
            return
        df = pd.read_csv(gwas_fp, sep='\t', usecols=cols, low_memory=False)
    df['SNP_ID_CURRENT'] = df['SNP_ID_CURRENT'].fillna('').apply(lambda x: clean_snp_ids(x))
    df = df.assign(SNPS=df['SNPS'].str.split(';')).explode('SNPS')
    df = df.assign(SNPS=df['SNPS'].str.split(',')).explode('SNPS')
    df['SNPS'] = df['SNPS'].str.strip()
    return df

def clean_snp_ids(snp_id):
    try:
        return f"rs{int(snp_id)}"
    except:
        if snp_id == '':
            return snp_id
        else:
            snp_id = snp_id.split('-')[0]
            snp_id = snp_id.split('_')[0]
            return f"rs{snp_id}" if not snp_id[:-1].isdigit() else f"rs{snp_id[:-1]}"

def find_disease(gwas, ppin_dir, out, ld, corr_thresh, window, population, ld_dir,
                 disable_pg=True, bootstrap=False):
    #print('Identifying GWAS traits.')
    sig_res = []
    eqtl_fps = [fp for fp in os.listdir(ppin_dir) if fp.endswith('snp_gene.txt')]
    # Total GWAS SNPs.
    M = gwas['SNPS'].nunique()
    for level_fp in eqtl_fps:
        level = f"{os.path.basename(level_fp).split('.')[0].split('_')[0]}"
        #print(f"\t{level}")
        df = pd.read_csv(os.path.join(ppin_dir, level_fp), sep='\t')
        snps = df['snp'].drop_duplicates().tolist()
        if ld:
            ld_snps = ld_proxy.ld_proxy(snps, corr_thresh, window, population, ld_dir, bootstrap)
            snps = ld_snps['rsidt'].drop_duplicates().tolist()
            df = (df.merge(ld_snps, left_on='snp', right_on='rsidq')
                  .sort_values(by=['snp', 'gene', 'dprime'])
                  .drop_duplicates(subset=['snp', 'gene', 'dprime']))
            write_results(df,  f'{level}_snp_gene.txt', out)
        overlap = gwas.loc[gwas['SNPS'].isin(snps)].drop_duplicates()
        # Total level eQTLs that are in the GWAS Catalog. 
        N = overlap['SNPS'].nunique()
        probs_df = []
        snp_trait_df = []
        for trait in tqdm(overlap['DISEASE/TRAIT'].drop_duplicates(),
                          disable=disable_pg):
            #print(f'\t\t{trait}')
            # Total trait-associated SNPs in GWAS Catalog 
            n = gwas[gwas['DISEASE/TRAIT']==trait]['SNPS'].nunique()
            trait_overlap = overlap[overlap['DISEASE/TRAIT'] == trait][['SNPS', 'DISEASE/TRAIT']]
            # Total trait-associated eQTLs
            X = trait_overlap['SNPS'].nunique()
            pval = hypergeom.sf(X-1, M, n, N)
            probs_df.append((level, trait, M, n, N, X, pval))
            snp_trait_df.append(trait_overlap.drop_duplicates())
        probs_cols = ['level', 'trait', 'total_gwas_snps', 'trait_snps',
                      'eqtls_in_catalog', 'trait_eqtls', 'pval']
        probs_df = pd.DataFrame(probs_df, columns=probs_cols)
        try:
            probs_df['adj_pval'] = mt.multipletests(probs_df['pval'], method='bonferroni')[1]
        except:
            continue
        sig_df = probs_df[probs_df['adj_pval'] < 0.05]
        sig_res.append(sig_df)
        snp_trait_df = pd.concat(snp_trait_df)
        snp_trait_df = snp_trait_df.rename(columns={'DISEASE/TRAIT': 'trait',
                                                    'SNPS': 'snp'})
        snp_trait_df = snp_trait_df[snp_trait_df['trait'].isin(sig_df['trait'])]
        if ld:
            write_results(
                snp_trait_df.merge(df.drop(columns=['snp']),
                                   how='inner', left_on='snp', right_on='rsidt'),
                f'{level}_sig_trait_snp_gene.txt', out)
        else:
            write_results(snp_trait_df.merge(df, how='inner'), f'{level}_sig_trait_snp_gene.txt', out)
        write_results(probs_df,  f'{level}_enrichment.txt', out)
    if len(sig_res) == 0:
        return
    sig_res = pd.concat(sig_res)
    write_results(sig_res, 'significant_enrichment.txt',  out)
    return sig_res

def write_results(res, fp, out):
    #print(f'\tWriting {fp}...')
    os.makedirs(out, exist_ok=True)
    res.to_csv(os.path.join(out, fp), sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Find disease associated with PPIN eQTLs.')
    parser.add_argument(
         '--ppin-dir', required=True,
        help='Filepath to directory containing eQTL-gene pairs for PPIN levels.')
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
        '--gwas', default=None,
        help='''Filepath to GWAS associations. 
        Default: Associations from the GWAS Catalog 
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '--ld', action='store_true', default=False,
        help='Include LD SNPs in identifying eQTLs and GWAS traits.')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched.')
    parser.add_argument('-p', '--population', default='EUR',
                        choices=['EUR'],
                        help='The ancestral population in which the LD is calculated.')
    parser.add_argument('--ld-dir', default='data/ld/dbs/super_pop/',
                        help='Directory containing LD database.')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    gwas = parse_gwas(args.gwas)
    find_disease(gwas, args.ppin_dir, args.output_dir,
                 args.ld, args.correlation_threshold, args.window, args.population, args.ld_dir)
    #write_results(res, args.output_dir)
    print('Done.')
