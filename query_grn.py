#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm

import non_spatial_eqtls
import ld_proxy

def parse_gwas(gwas_fp):
    print('Parsing GWAS associations...')
    cols = ['SNPS', 'SNP_ID_CURRENT', 'DISEASE/TRAIT', 'PUBMEDID',
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

def parse_snp_input(snp_fp):
    if os.path.isfile(snp_fp[0]):
        df = pd.read_csv(snp_fp[0], sep='\t', header=None, names=['snp'])
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

def extract_trait_snps(trait, gwas):
    return gwas[gwas['DISEASE/TRAIT'] == trait]['SNPS'].drop_duplicates()

def extract_pmid_snps(pmid, gwas):
    return gwas[gwas['PUBMEDID'] == int(pmid)]['SNPS'].drop_duplicates()

def parse_grn(grn_dir):
    print('Parsing tissue gene regulatory map...')
    grn = []
    for chrom in os.listdir(grn_dir):
        chrom_eqtls = pd.read_csv(
            os.path.join(grn_dir, chrom, 'significant_eqtls.txt'), sep='\t')
        grn.append(chrom_eqtls)
    return pd.concat(grn)
    
def get_eqtls(snps, grn, output_dir, non_spatial, non_spatial_dir, snp_ref_dir, gene_ref_dir,
ld, corr_thresh, window, population, ld_dir):
    if ld:
        ld_snps = ld_proxy.ld_proxy(snps, corr_thresh, window, population, ld_dir)
        snps = ld_snps['rsidt'].drop_duplicates()
        write_results(ld_snps, os.path.join(output_dir, 'query_snp_ld.txt'))
    constrained_eqtls = grn[grn['snp'].isin(snps)]
    cols = cols = ['snp', 'gene', 'gencode_id', 'tissue', 'interaction_type',
                   'eqtl_type','gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    if non_spatial:
        tissue = grn['tissue'].drop_duplicates()[0]
        unconstrained_eqtls = non_spatial_eqtls.get_eqtls(snps, tissue, output_dir, non_spatial_dir,
                                                        snp_ref_dir, gene_ref_dir)
    if len(constrained_eqtls) > 0:
        constrained_eqtls['eqtl_type'] = 'spatial'
    eqtls = pd.concat([constrained_eqtls, unconstrained_eqtls])
    if len(eqtls) == 0:
        sys.exit('No eQTLs found in the gene regulatory map.')
    write_results(eqtls, os.path.join(output_dir, 'query_eqtls.txt'))
    return eqtls

def write_results(res, fp):
    output_dir = os.path.dirname(fp)
    os.makedirs(output_dir, exist_ok=True)
    res.to_csv(fp, sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Find disease associated with PPIN eQTLs.')
    parser.add_argument(
        '-s', '--snps',
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
        '--non-spatial-dir', default='data/GTEx/', help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '-g', '--gwas', default=None,
        help='''Filepath to GWAS associations. 
        Default: Associations from the GWAS Catalog 
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '--snp-ref-dir', default='data/snps/', help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default='data/genes/', help='Filepath to gene BED.')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    if not args.snps and not args.trait and not args.pmid:
        sys.exit('FATAL: One of --snps, --trait, or --pmid is required.\nExiting.')
    gwas = parse_gwas(args.gwas)
    snps = args.snps
    print(snps)
    if args.trait:
        snps = extract_trait_snps(args.trait, gwas)
    elif args.pmid:
        snps = extract_pmid_snps(args.pmid, gwas)
    else:
        snps = parse_snp_input(args.snps)
    grn = parse_grn(args.grn_dir)
    eqtls = get_eqtls(snps, grn, args.output_dir, args.non_spatial, args.non_spatial_dir,
                      args.snp_ref_dir, args.gene_ref_dir)
    
