#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm
from time import time
import multiprocessing as mp
from itertools import repeat
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool

    
def get_eqtls(snps, tissue, output_dir, non_spatial_dir, snp_ref_dir, gene_ref_dir):
    snp_df, failed = rsids2pos(snps, snp_ref_dir)
    if len(failed) > 0:
        merged_rsids = merged_snps(failed, snp_ref_dir)
        merged_df, failed = rsids2pos(merged_rsids['old'].tolist(), snp_ref_dir)
        snp_df = pd.concat([snp_df, merged_df])
        if len(merged_rsids) > 0:
            num_merged_rsids = merged_rsids['new'].nunique()
            write_file(merged_rsids, os.path.join(output_dir, 'merged_snps.txt'))
            print(f'WARNING: {num_merged_rsids} SNPs have been merged. '+
                  'See "merged_snps.txt" for details')
        if len(failed) > 0:
            write_file(pd.DataFrame({'failed_snps': failed}),
                       os.path.join(output_dir, 'failed_snps.txt'))
            print(f'''WARNING: {len(failed)} SNPs could not be found in the SNP database.''')
    snp_df['snp_locus'] = snp_df['start'] + 1
    snp_df['id'] = 'chr' + snp_df['chrom'].astype(str) + '_' + snp_df['snp_locus'].astype(str)
    snp_df = snp_df.rename(columns={'chrom': 'snp_chr',
                                    'rsid': 'snp'})
    cis_eqtls = get_cis_eqtls(snp_df, tissue, non_spatial_dir, gene_ref_dir)
    trans_eqtls = get_trans_eqtls(snp_df, tissue, non_spatial_dir, gene_ref_dir)
    eqtls = pd.concat([cis_eqtls, trans_eqtls])
    if len(eqtls) == 0:
        sys.exit('No eQTLs found in the gene regulatory map.')
    write_file(eqtls, os.path.join(output_dir, 'non_spatial_eqtls.txt'))        
    return eqtls

def get_trans_eqtls(snps, tissue, non_spatial_dir, gene_ref_dir):
    fp = os.path.join(non_spatial_dir, 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')
    chunk_size = 200000
    res = []
    cols = ['variant_id', 'snp', 'gencode_id', 'gene', 'tissue', 'interaction_type', 'eqtl_type',
            'tss_distance', 'maf', 'pval_nominal','slope', 'slope_se',
            'gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    for df in pd.read_csv(fp, sep='\t', chunksize=chunk_size):
        df['id'] = df['variant_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
        df = df.rename(columns = {'tissue_id': 'tissue',
                                  'tissue_af': 'maf',
                                  'gene_name': 'gene'})
        res.append(df[((df['id'].isin(snps['id'])) & (df['tissue'] == tissue))])
    if len(res) == 0:
        return
    res = pd.concat(res)
    res = res.rename(columns={'gene_id': 'gencode_id'})
    res = res.merge(snps, how='left', on='id')
    gene_ref = pd.read_csv(f'{gene_ref_dir}gene_reference.bed', sep='\t', header=None,
                           names=['gene_chr', 'gene_start', 'gene_end', 'gene', 'gencode_id'])
    res = res.merge(gene_ref, how='left', on=['gene','gencode_id', 'gene_chr'])
    res['interaction_type'] = 'Trans-intrachromosomal'
    res.loc[res['gene_chr'] != res['snp_chr'], 'interaction_type'] = 'Trans-interchromosomal'
    res['eqtl_type'] = 'non_spatial'
    res['tss_distance'] = -1
    return res[cols]

def get_cis_eqtls(snps, tissue, non_spatial_dir, gene_ref_dir):
    cis_dir = os.path.join(non_spatial_dir, 'GTEx_Analysis_v8_eQTL')
    tissue_fp = [fp for fp in os.listdir(cis_dir) 
               if fp.endswith('gene_pairs.txt.gz') and
               fp.split('.')[0] == tissue][0]
    chunk_size = 200000
    res = []
    cols = ['variant_id', 'snp', 'gencode_id', 'gene', 'tissue', 'interaction_type', 'eqtl_type',
            'tss_distance', 'maf', 'pval_nominal','slope', 'slope_se',
            'gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    for df in pd.read_csv(f'{cis_dir}/{tissue_fp}', sep='\t',
                          chunksize=chunk_size, compression='gzip'):
        df['id'] = df['variant_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
        res.append(df[df['id'].isin(snps['id'])])
    if len(res) == 0:
        return
    res = pd.concat(res)
    res = res.rename(columns={'gene_id': 'gencode_id'})
    res = res.merge(snps, how='left', on='id')
    gene_ref = pd.read_csv(f'{gene_ref_dir}gene_reference.bed', sep='\t', header=None,
                           names=['gene_chr', 'gene_start', 'gene_end', 'gene', 'gencode_id'])
    res = res.merge(gene_ref, how='left', on='gencode_id')
    res['tissue'] = tissue
    res['interaction_type'] = 'Cis'
    res['eqtl_type'] = 'non_spatial'
    return res[cols]
    
def write_file(df, fp):
    dirname = os.path.dirname(fp)
    os.makedirs(dirname, exist_ok=True)
    df.to_csv(fp, sep='\t', index=False)

def parse_snp_input(snp_fp):
    if os.path.isfile(snp_fp[0]):
        df = pd.read_csv(snp_fp[0], sep='\t', header=None, names=['snp'])
        return df[df['snp'] != "snp"]['snp'].drop_duplicates().tolist()
    else:
        return list(set(snp_fp))

def rsids2pos(snps, ref_dir):
    chrom_list = [fp.split('.')[0].split('_')[-1]
                  for fp in os.listdir(f'{ref_dir}/dbs/')
                  if fp.startswith('bed_chr')]
    manager = mp.Manager()
    chrom_dict = manager.dict()
    for chrom in chrom_list:
        chrom_dict[chrom] = pd.DataFrame()
    '''
    for chrom in chrom_list:
        rsids2pos_chrom(chrom, snps, chrom_dict, ref_dir)
    '''
    with mp.Pool(16) as pool:
        pool.starmap(rsids2pos_chrom,
                     zip(chrom_list,
                        repeat(snps),
                        repeat(chrom_dict),
                        repeat(ref_dir)))
    df = []
    for chrom in chrom_dict.keys():
        df.append(chrom_dict[chrom])
    df = pd.concat(df)
    missed = list(set(snps).difference(set(df['rsid'].tolist())))
    
    return df, missed

def merged_snps(query_snps, ref_dir):
    pd.options.mode.chained_assignment = None
    merged_archive = 'human_9606_b151_GRCh38p7_RsMergeArch.bcp.db'
    db = create_engine(f'sqlite:///{ref_dir}/dbs/{merged_archive}', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM merged_arch WHERE new = {}; '''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            res = pd.read_sql(sql.format(snp[2:]), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
    snps['new'] = 'rs' + snps['new'].astype(str)
    snps['old'] = 'rs' + snps['old'].astype(str)
    return snps
    

def rsids2pos_chrom(chrom, query_snps, chrom_dict, ref_dir):
    pd.options.mode.chained_assignment = None
    db = create_engine(f'sqlite:///{ref_dir}/dbs/bed_chr_{chrom}.db', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM snps WHERE rsid = '{}'; '''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            res = pd.read_sql(sql.format(snp), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
    chrom_dict[chrom] = snps

def parse_args():
    parser = argparse.ArgumentParser(
        description='Find disease associated with PPIN eQTLs.')
    parser.add_argument(
        '-s', '--snps', nargs='+',
        help='''A space-separated list of SNP rsIDs or filepath to a file 
        containing query SNP rsIDs in the 'snp' column. Note: this flag is mutually 
        exclusive with the --trait and --pmid flag.''')
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
        '-t', '--tissue', required=True, help='Tissue in which eQTLs are mapped.')
    parser.add_argument(
        '--non-spatial-dir', default='data/GTEx/', help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '--snp-ref-dir', default='data/snps/', help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default='data/genes/', help='Filepath to gene BED.')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    snps = parse_snp_input(args.snps)
    eqtls = get_eqtls(snps, args.tissue, args.output_dir,
                      args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir)










