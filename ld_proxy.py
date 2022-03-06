#!/usr/bin/env python
import pandas as pd
from time import time
import multiprocessing as mp
from itertools import repeat
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import sys
import os
import argparse

def parse_snps(snps_input):
    if os.path.isfile(snps_input[0]): 
        df = pd.read_csv(snps_input[0], sep='\t', header=None, names=['snp'])
        df = df[df['snp'] != 'snp']
        return df.drop_duplicates()
    else:
        return pd.DataFrame({'snp': snps_input})

def ld_proxy_chrom(chrom, query_snps, chrom_dict, corr_thresh, window, pop, ld_dir):
    pd.options.mode.chained_assignment = None
    db = create_engine(f'sqlite:///{ld_dir}/{chrom}.db', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM ld WHERE rsidq = '{}'; '''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            res = pd.read_sql(sql.format(snp), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
    snps = snps[((snps['corr'] >= corr_thresh) & 
                              (abs(snps['posq']-snps['post']) <= window))]
    snps.loc[:, 'chromt_post'] = snps['chromt'].astype(str)  + '_' + snps['post'].astype(str) 
    res = snps[['rsidt', 'chromt_post']]
    res = pd.DataFrame(res['rsidt'].str.split(';').tolist(), index=res['chromt_post']).stack()
    res = res.reset_index()[[0, 'chromt_post']] 
    res.columns = ['rsidt_res', 'chromt_post']
    snps = (snps
                .merge(res, how='inner', on='chromt_post')
               .drop(columns=['chromt_post', 'rsidt'])
               .rename(columns={'rsidt_res': 'rsidt'})
           )
    chrom_dict[chrom] = snps[['chromq',	'posq', 'rsidq', 'chromt', 'post', 'rsidt', 'corr','dprime']]

def ld_proxy(query_snps, corr_thresh, window, pop, ld_dir):
    ld_dir = os.path.join(ld_dir, pop)

    chrom_list = [fp.split('.')[0] for fp in os.listdir(ld_dir) if fp.startswith('chr')]
    manager = mp.Manager()
    chrom_dict = manager.dict()
    for chrom in chrom_list:
        chrom_dict[chrom] = pd.DataFrame()
    with mp.Pool(16) as pool:
        pool.starmap(ld_proxy_chrom,
                     zip(chrom_list,
                        repeat(query_snps),
                        repeat(chrom_dict),
                        repeat(corr_thresh),
                        repeat(window),
                        repeat(pop),
                        repeat(ld_dir)))
    df = []
    for chrom in chrom_dict.keys():
        df.append(chrom_dict[chrom])
    df = pd.concat(df)
    query_snps_df = pd.DataFrame({'chromq': '',
                        'posq': '',
                        'rsidq': query_snps,
                        'chromt': '',
                        'post': '',
                        'rsidt': query_snps,
                        'corr': 1,
                        'dprime': 1
                       })
    df = (pd.concat([df, query_snps_df])
             .sort_values(by=['rsidq', 'corr', 'dprime'], ascending=False))
    return df

def write_results(df, out):
    #print('Writing PPIN...')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep='\t', index=False)
    
def parse_args():
    parser = argparse.ArgumentParser(
        description='A tool to find SNPs in LD.')
    parser.add_argument(
        '-s', '--snps', required=True, nargs='+',
        help='''A space-separated list of rsIDs or filepath to a file 
        containing SNP rsIDs in one column.''' )
    parser.add_argument(
        '-o', '--output', required=True, help='Filepath to write results.')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched.')
    parser.add_argument('-p', '--population', default='EUR',
                        choices=['EUR'],
                        help='The ancestral population in which the LD is calculated')
    parser.add_argument('--ld-dir', default='data/ld/dbs/super_pop/',
                        help='Directory containing LD database.')
    return parser.parse_args()


if __name__=='__main__':
    start_time = time()
    args = parse_args()
    snps = parse_snps(args.snps)
    df = ld_proxy(snps['snp'].tolist(), args.correlation_threshold,
                  args.window, args.population, args.ld_dir)
    write_results(df, args.output)
    print('Total time elasped: {:.2f} mins.'.format(
        (time()-start_time)/60))

