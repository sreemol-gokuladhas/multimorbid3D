#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm
import time
import ssl

import non_spatial_eqtls
import ld_proxy
import logger

def parse_gwas(gwas_fp, logger):
    logger.write('Parsing GWAS associations...')
    cols = ['SNPS', 'SNP_ID_CURRENT', 'DISEASE/TRAIT', 'PUBMEDID',
            'P-VALUE',  'OR or BETA', '95% CI (TEXT)']
    if gwas_fp is None:
        fp = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'
        # To prevent the occasssional urllib.error.URLError
        ssl._create_default_https_context = ssl._create_unverified_context
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

def parse_snp_input(snp_fp, logger):
    if os.path.isfile(snp_fp[0]):
        df = pd.read_csv(snp_fp[0], sep='\t', header=None, names=['snp'])
        df = df[df['snp'] != "snp"]
        df['snp'] = df['snp'].str.strip()
        return df[df['snp'] != "snp"]['snp'].drop_duplicates().tolist()
    else:
        return list(set(snp_fp))
    
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

def extract_trait_snps(trait, gwas, logger):
    logger.write(f'Extracting SNPs from trait entries with the word "{trait}"')
    return gwas[gwas['DISEASE/TRAIT'] == trait]['SNPS'].drop_duplicates()

def extract_pmid_snps(pmid, gwas, logger):
    logger.write(f'Extracting SNPs from PubMed ID {pmid}')
    return gwas[gwas['PUBMEDID'] == int(pmid)]['SNPS'].drop_duplicates()

def parse_grn(grn_dir, logger):
    logger.write('Parsing tissue gene regulatory map...')
    grn = []
    chrom_dirs =  [d for d in os.listdir(grn_dir) if os.path.isdir(os.path.join(grn_dir, d))]#not d.endswith('.log')]
    for chrom in chrom_dirs:
        chrom_eqtls = pd.read_csv(
            os.path.join(grn_dir, chrom, 'significant_eqtls.txt'), sep='\t')
        grn.append(chrom_eqtls)
    return pd.concat(grn)
    
def get_eqtls(snps, grn, output_dir,
              non_spatial, non_spatial_dir, snp_ref_dir, gene_ref_dir,
              ld, corr_thresh, window, population, ld_dir, logger, bootstrap=False):
    if ld:
        ld_snps = ld_proxy.ld_proxy(snps, corr_thresh, window, population, ld_dir, logger, bootstrap)
        snps = ld_snps['rsidt'].drop_duplicates()
        write_results(ld_snps, os.path.join(output_dir, 'query_snp_ld.txt'))
    constrained_eqtls = grn[grn['snp'].isin(snps)]
    unconstrained_eqtls = pd.DataFrame()
    cols = cols = ['snp', 'gene', 'gencode_id', 'tissue', 'interaction_type',
                   'eqtl_type','gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    if non_spatial:
        tissue = grn['tissue'].drop_duplicates()[0]
        unconstrained_eqtls = non_spatial_eqtls.get_eqtls(
            snps, tissue, output_dir, non_spatial_dir,
            snp_ref_dir, gene_ref_dir, logger)
    if len(constrained_eqtls) > 0:
        constrained_eqtls.loc[:, 'eqtl_type'] = 'spatial'
    eqtls = pd.concat([constrained_eqtls, unconstrained_eqtls])
    if eqtls.empty:
        logger.write('No eQTLs found in the gene regulatory map.')
        return pd.DataFrame()
    write_results(eqtls.drop_duplicates(), os.path.join(output_dir, 'query_eqtls.txt'))
    return eqtls

def get_gene_eqtls(gene_list, grn, output_dir,
                   non_spatial, non_spatial_dir, snp_ref_dir, gene_ref_dir,
                   logger, bootstrap=False):
    #logger.write('Identifying gene eQTLs...')
    res = []
    for level in range(len(gene_list)):
        constrained_eqtls = (
            grn[grn['gene'].isin(gene_list[level])][['snp', 'gene']]
            .drop_duplicates()
            .assign(eqtl_type = 'spatial'))
        unconstrained_eqtls = pd.DataFrame()
        if non_spatial:
            tissue = grn['tissue'].drop_duplicates()[0]
            unconstrained_eqtls = non_spatial_eqtls.get_gene_eqtls(
                gene_list[level], tissue, output_dir,
                non_spatial_dir, snp_ref_dir, gene_ref_dir, logger, bootstrap=bootstrap)
            if not unconstrained_eqtls.empty:
                unconstrained_eqtls = (unconstrained_eqtls[['snp', 'gene']]
                                       .drop_duplicates()
                                       .assign(eqtl_type = 'non_spatial'))
        df = pd.concat([constrained_eqtls, unconstrained_eqtls])
        df.loc[df.duplicated(subset=['snp', 'gene'], keep=False), 'eqtl_type'] = 'both'
        df = df.drop_duplicates()
        write_results(df,
                      os.path.join(output_dir, f'level{level}_snp_gene.txt'))
        del df
        #res.append(df)
    #return pd.concat(res)
        

def write_results(res, fp):
    output_dir = os.path.dirname(fp)
    os.makedirs(output_dir, exist_ok=True)
    res.to_csv(fp, sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Query tissue GRN for eQTL associations.')
    parser.add_argument(
        '-s', '--snps', nargs='+',
        help='''A space-separated list of SNP rsIDs or filepath to a file 
        containing query SNP rsIDs in the 'snp' column. Note: this flag is mutually 
        exclusive with the --trait and --pmid flag.''')
    parser.add_argument(
        '-t', '--trait',
        help='''GWAS trait to query. Note: this flag is mutually exclusive with
        the --snps and --pmid flags''')
    parser.add_argument(
        '-p', '--pmid',
        help='''PubMed ID of the GWAS to query. Note: this flag is mutually exclusive with
        the --snps and --trait flag''')
    parser.add_argument(
        '--grn-dir', required=True,
        help='''Directory containing tissue gene regulatory network. 
        The subdirectories should contain significant_eqtls.txt for each chromosome.''')
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
        '--non-spatial', action='store_true', default=False,
        help='Include non-spatial eQTLs.')
    parser.add_argument(
        '--non-spatial-dir', default=os.path.join(os.path.dirname(__file__), 'data/GTEx/'),
        help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '-g', '--gwas', default=None,
        help='''Filepath to GWAS associations. 
        Default: Associations from the GWAS Catalog 
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '--snp-ref-dir', default=os.path.join(os.path.dirname(__file__), 'data/snps/'),
        help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default=os.path.join(os.path.dirname(__file__),'data/genes/'),
        help='Filepath to gene BED.')
    parser.add_argument(
        '--ld', action='store_true', default=False,
        help='Include LD SNPs in identifying eQTLs and GWAS traits.')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched.')
    parser.add_argument(
        '--population', default='EUR', choices=['EUR'],
        help='The ancestral population in which the LD is calculated.')
    parser.add_argument(
        '--ld-dir',
        default=os.path.join(os.path.dirname(__file__), 'data/ld/dbs/super_pop/'),
        help='Directory containing LD database.')
    return parser.parse_args()


if __name__=='__main__':
    pd.options.mode.chained_assignment = None
    args = parse_args()
    if not args.snps and not args.trait and not args.pmid:
        sys.exit('FATAL: One of --snps, --trait, or --pmid is required.\nExiting.')
    start_time = time.time()
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'query_grn.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    gwas = parse_gwas(args.gwas, logger)
    snps = args.snps
    if args.trait:
        snps = extract_trait_snps(args.trait, gwas, logger)
    elif args.pmid:
        snps = extract_pmid_snps(args.pmid, gwas, logger)
    else:
        snps = parse_snp_input(args.snps, logger)
    grn = parse_grn(args.grn_dir, logger)
    eqtls = get_eqtls(snps, grn, args.output_dir, args.non_spatial, args.non_spatial_dir,
                      args.snp_ref_dir, args.gene_ref_dir,
                      args.ld, args.correlation_threshold, args.window, args.population,
                      args.ld_dir, logger)
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
