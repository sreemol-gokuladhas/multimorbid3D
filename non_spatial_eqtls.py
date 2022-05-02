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

import logger
    
def get_eqtls(snps, tissue, output_dir, non_spatial_dir, snp_ref_dir, gene_ref_dir,
              logger):
    snp_df, failed = rsids2pos(snps, snp_ref_dir)
    if len(failed) > 0:
        merged_rsids = merged_snps(failed, snp_ref_dir)
        if not merged_rsids.empty:
            merged_df, failed = rsids2pos(merged_rsids['old'].tolist(), snp_ref_dir)
            if not merged_df.empty:
                snp_df = pd.concat([snp_df, merged_df])
            if len(merged_rsids) > 0:
                num_merged_rsids = merged_rsids['new'].nunique()
                write_file(merged_rsids, os.path.join(output_dir, 'merged_snps.txt'))
                logger.write(f'WARNING: {num_merged_rsids} SNPs have been merged. '+
                             'See "merged_snps.txt" for details')
        if len(failed) > 0:
            write_file(pd.DataFrame({'failed_snps': failed}),
                       os.path.join(output_dir, 'failed_snps.txt'))
            logger.write(f'''WARNING: {len(failed)} SNPs could not be found in the SNP database.''')
    snp_df['snp_locus'] = snp_df['start'] + 1
    snp_df['id'] = 'chr' + snp_df['chrom'].astype(str) + '_' + snp_df['snp_locus'].astype(str)
    snp_df = snp_df.rename(columns={'chrom': 'snp_chr',
                                    'rsid': 'snp'})
    cis_eqtls = get_cis_eqtls(snp_df, tissue, non_spatial_dir, gene_ref_dir)
    trans_eqtls = get_trans_eqtls(snp_df, tissue, non_spatial_dir, gene_ref_dir)
    eqtls = pd.concat([cis_eqtls, trans_eqtls])
    if eqtls.empty:
        return pd.DataFrame()
    #sys.exit('No eQTLs found in the gene regulatory map.')
    #write_file(eqtls, os.path.join(output_dir, 'non_spatial_eqtls.txt'))        
    return eqtls

def get_trans_eqtls(snps, tissue, non_spatial_dir, gene_ref_dir):
    fp = os.path.join(non_spatial_dir, 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')
    chunk_size = 200000
    res = []
    cols = ['variant_id', 'snp', 'gencode_id', 'gene', 'tissue', 'interaction_type', 'eqtl_type',
            'tss_distance', 'maf', 'pval_nominal','slope', 'slope_se',
            'gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    for df in pd.read_csv(fp, sep='\t', chunksize=chunk_size):
        df['chr'] = df['variant_id'].apply(lambda x: x.split('_')[0]) 
        df['pos'] = df['variant_id'].apply(lambda x: int(x.split('_')[1]))
        df['ref'] = df['variant_id'].apply(lambda x: len(x.split('_')[2])) 
        df['alt'] = df['variant_id'].apply(lambda x: len(x.split('_')[3]))
        # Insertions have one base added.
        df.loc[df['ref'] > 1, 'pos'] = df['pos'].astype(int) + 1
        df['id'] = df['chr'] + '_' + df['pos'].astype(str)
        #df['id'] = df['variant_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
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
        df['chr'] = df['variant_id'].apply(lambda x: x.split('_')[0]) 
        df['pos'] = df['variant_id'].apply(lambda x: int(x.split('_')[1]))
        df['ref'] = df['variant_id'].apply(lambda x: len(x.split('_')[2])) 
        df['alt'] = df['variant_id'].apply(lambda x: len(x.split('_')[3]))
        # Insertions have one base added.
        df.loc[df['ref'] > 1, 'pos'] = df['pos'].astype(int) + 1
        df['id'] = df['chr'] + '_' + df['pos'].astype(str)
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


def get_gene_eqtls(genes, tissue, output_dir, non_spatial_dir, snp_ref_dir, gene_ref_dir,
                   logger, bootstrap=False):
    cis_eqtls = get_cis_gene_eqtls(genes, tissue, non_spatial_dir, gene_ref_dir,
                                   snp_ref_dir, bootstrap)
    trans_eqtls = get_trans_gene_eqtls(genes, tissue, non_spatial_dir, gene_ref_dir,
                                       snp_ref_dir, bootstrap)
    eqtls = pd.concat([cis_eqtls, trans_eqtls])
    if eqtls.empty:
        #logger.write('No eQTLs found in the gene regulatory map.')
        return pd.DataFrame()
    #write_file(eqtls, os.path.join(output_dir, 'non_spatial_eqtls.txt'))
    return eqtls


def get_trans_gene_eqtls(genes, tissue, non_spatial_dir, gene_ref_dir, snp_ref_dir,
                         bootstrap):
    fp = os.path.join(non_spatial_dir, 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')
    chunk_size = 200000
    res = []
    cols = ['variant_id', 'snp', 'gencode_id', 'gene', 'tissue', 'interaction_type', 'eqtl_type',
            'tss_distance', 'maf', 'pval_nominal','slope', 'slope_se',
            'gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    gene_ref = pd.read_csv(f'{gene_ref_dir}gene_reference.bed', sep='\t', header=None,
                           names=['gene_chr', 'gene_start', 'gene_end', 'gene', 'gencode_id'])
    gene_ref['id'] = gene_ref['gencode_id'].str[:15]
    for df in pd.read_csv(fp, sep='\t', chunksize=chunk_size):
        df = df.rename(columns = {'tissue_id': 'tissue',
                                  'tissue_af': 'maf',
                                  'gene_name': 'gene'})
        #df['id'] = df['gene_id'].str[:15]
        df = df.merge(gene_ref, how='inner', on=['gene', 'gene_chr'])
        res.append(df[((df['gene'].isin(genes)) & (df['tissue'] == tissue))])
    if len(res) == 0:
        return pd.DataFrame()
    res = pd.concat(res)
    if res.empty:
        return pd.DataFrame()
    res['snp_chr'] = res['variant_id'].apply(lambda x: x.split('_')[0])
    res['snp_locus'] = res['variant_id'].apply(lambda x: x.split('_')[1])
    snps = pos2rsids(res[['snp_chr', 'snp_locus']].drop_duplicates(),
                     snp_ref_dir, bootstrap)
    res['snp_locus'] = res['snp_locus'].astype(int)
    #snp['snp_locus'] = snp['snp_locus'].astype(int)
    res = res.merge(snps, how='inner', on=['snp_chr', 'snp_locus'])
    res['tissue'] = tissue
    res['interaction_type'] = 'Cis'
    res['eqtl_type'] = 'non_spatial'
    res['interaction_type'] = 'Trans-intrachromosomal'
    res.loc[res['gene_chr'] != res['snp_chr'], 'interaction_type'] = 'Trans-interchromosomal'
    res['eqtl_type'] = 'non_spatial'
    res['tss_distance'] = -1
    return res[cols]

def get_cis_gene_eqtls(genes, tissue, non_spatial_dir, gene_ref_dir,
                       snp_ref_dir, bootstrap):
    cis_dir = os.path.join(non_spatial_dir, 'GTEx_Analysis_v8_eQTL')
    tissue_fp = [fp for fp in os.listdir(cis_dir) 
               if fp.endswith('gene_pairs.txt.gz') and
               fp.split('.')[0] == tissue][0]
    chunk_size = 200000
    res = []
    cols = ['variant_id', 'snp', 'gencode_id', 'gene', 'tissue', 'interaction_type', 'eqtl_type',
            'tss_distance', 'maf', 'pval_nominal','slope', 'slope_se',
            'gene_chr', 'gene_start', 'gene_end', 'snp_chr', 'snp_locus']
    gene_ref = pd.read_csv(f'{gene_ref_dir}gene_reference.bed', sep='\t', header=None,
                           names=['gene_chr', 'gene_start', 'gene_end', 'gene', 'gencode_id'])
    gene_ref['id'] = gene_ref['gencode_id'].str[:15]
    for df in pd.read_csv(f'{cis_dir}/{tissue_fp}', sep='\t',
                          chunksize=chunk_size, compression='gzip'):
        df['id'] = df['gene_id'].str[:15]
        df = df.merge(gene_ref, how='inner', on='id')
        res.append(df[df['gene'].isin(genes)])

    if len(res) == 0:
        return pd.DataFrame()
    res = pd.concat(res)
    if res.empty:
        return pd.DataFrame()
    res['snp_chr'] = res['variant_id'].apply(lambda x: x.split('_')[0])
    res['snp_locus'] = res['variant_id'].apply(lambda x: int(x.split('_')[1]))
    snps = pos2rsids(res[['snp_chr', 'snp_locus']].drop_duplicates(),
                     snp_ref_dir, bootstrap)
    res = res.merge(snps, how='inner', on=['snp_chr', 'snp_locus'])
    res['tissue'] = tissue
    res['interaction_type'] = 'Cis'
    res['eqtl_type'] = 'non_spatial'
    return res[cols]
    
def write_file(df, fp):
    dirname = os.path.dirname(fp)
    os.makedirs(dirname, exist_ok=True)
    df.to_csv(fp, sep='\t', index=False)

def parse_snp_input(snp_fp, logger):
    if os.path.isfile(snp_fp[0]):
        df = pd.read_csv(snp_fp[0], sep='\t', header=None, names=['snp'])
        return df[df['snp'] != "snp"]['snp'].drop_duplicates().tolist()
    else:
        return list(set(snp_fp))

def rsids2pos(snps, ref_dir, bootstrap=False):
    chrom_list = [fp.split('.')[0].split('_')[-1]
                  for fp in os.listdir(f'{ref_dir}/dbs/')
                  if fp.startswith('bed_chr')]
    manager = mp.Manager()
    chrom_dict = manager.dict()
    for chrom in chrom_list:
        chrom_dict[chrom] = pd.DataFrame()
    if bootstrap:
        for chrom in chrom_list:
            rsids2pos_chrom(chrom, snps, chrom_dict, ref_dir)
    else:
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


def rsids2pos_chrom(chrom, query_snps, chrom_dict, ref_dir):
    pd.options.mode.chained_assignment = None
    db = create_engine(f'sqlite:///{ref_dir}/dbs/bed_chr_{chrom}.db', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM snps WHERE rsid = '{}'; '''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            res = pd.read_sql(sql.format(snp.strip()), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) > 0:
        snps = pd.concat(snps)
        chrom_dict[chrom] = snps

def merged_snps(query_snps, ref_dir):
    pd.options.mode.chained_assignment = None
    merged_archive = 'human_9606_b151_GRCh38p7_RsMergeArch.bcp.db'
    db = create_engine(f'sqlite:///{ref_dir}/dbs/{merged_archive}', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM merged_arch WHERE new = {}; '''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            try:
                res = pd.read_sql(sql.format(snp[2:]), conn)
                if not res.empty:
                    snps.append(res)
            except:
                pass
    if len(snps) == 0:
        return pd.DataFrame()
    snps = pd.concat(snps)
    snps['new'] = 'rs' + snps['new'].astype(str)
    snps['old'] = 'rs' + snps['old'].astype(str)
    return snps


def pos2rsids(snps, ref_dir, bootstrap=False):
    chrom_list = [fp.split('.')[0].split('_')[-1]
                  for fp in os.listdir(f'{ref_dir}/dbs/')
                  if fp.startswith('bed_chr')]
    if bootstrap:
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom] = pd.DataFrame()
        for chrom in chrom_list:
            pos2rsids_chrom(chrom, snps, chrom_dict, ref_dir)
    else:
        manager = mp.Manager()
        chrom_dict = manager.dict()
        for chrom in chrom_list:
            chrom_dict[chrom] = pd.DataFrame()
        with mp.Pool(16) as pool:
            pool.starmap(pos2rsids_chrom,
                         zip(chrom_list,
                            repeat(snps),
                            repeat(chrom_dict),
                            repeat(ref_dir)))
    df = []
    for chrom in chrom_dict.keys():
        df.append(chrom_dict[chrom])
    df = pd.concat(df)
    if df.empty:
        return pd.DataFrame()
    df = df.rename(columns={'chrom': 'snp_chr',
                            'locus': 'snp_locus',
                            'rsid': 'snp'})
    df['snp_chr'] = 'chr' + df['snp_chr'].astype(str)
    return df


def pos2rsids_chrom(chrom, query_snps, chrom_dict, ref_dir):
    pd.options.mode.chained_assignment = None
    db = create_engine(f'sqlite:///{ref_dir}/dbs/bed_chr_{chrom}.db', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM snps WHERE chrom = '{}' and start = {}; '''
    snps = []
    query = query_snps[query_snps['snp_chr'].str[3:] == chrom]
    if query.empty:
        return
    with db.connect() as conn:
        for _, snp in query.iterrows():
            #logger.write(snp)
            res = pd.read_sql(sql.format(snp['snp_chr'][3:],
                                         int(snp['snp_locus']) - 1), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
    snps['locus'] = snps['start'] + 1
    chrom_dict[chrom] = snps

def parse_args():
    parser = argparse.ArgumentParser(
        description='Identify non-spatial eQTLs.')
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
        '--non-spatial-dir', default=os.path.join(os.path.dirname(__file__), 'data/GTEx/'),
        help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '--snp-ref-dir', default=os.path.join(os.path.dirname(__file__), 'data/snps/'),
        help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default=os.path.join(os.path.dirname(__file__), 'data/genes/'),
        help='Filepath to gene BED.')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'non_spatial_eqtls.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    snps = parse_snp_input(args.snps, logger)
    eqtls = get_eqtls(snps, args.tissue, args.output_dir,
                      args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir,
                      logger)
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')










