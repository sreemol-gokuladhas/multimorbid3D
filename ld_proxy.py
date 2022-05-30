#!/usr/bin/env python
import pandas as pd
from time import time
import multiprocessing as mp
from itertools import repeat
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import sys
import os
import argparse
from tqdm import tqdm

import logger
import non_spatial_eqtls as ns

def parse_snps(snps_input, logger):
    if os.path.isfile(snps_input[0]): 
        df = pd.read_csv(snps_input[0], sep='\t', header=None, names=['snp'])
        df = df[~df['snp'].str.lower().isin(['snp', 'snps'])]
        return df.drop_duplicates()
    else:
        return pd.DataFrame({'snp': snps_input})
    
    
def ld_proxy_chrom(chrom, query_snps, corr_thresh, window, pop, ld_dir):
    pd.options.mode.chained_assignment = None
    db = create_engine(f'sqlite:///{ld_dir}/{chrom}.db', 
                       echo=False, poolclass=NullPool)
    sql = '''SELECT * FROM ld WHERE rsidq = '{}' AND corr >= {};'''
    snps = []
    with db.connect() as conn:
        for snp in query_snps:
            res = pd.read_sql(sql.format(snp, corr_thresh), conn)
            if not res.empty:
                snps.append(res)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
    snps = snps[(abs(snps['posq']-snps['post']) <= window)]
    if snps.empty:
        return
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



def ld_proxy(query_snps, corr_thresh, window, pop, ld_dir, logger, bootstrap):
    ld_dir = os.path.join(ld_dir, pop)
    chrom_list = [fp.split('.')[0] for fp in os.listdir(ld_dir) if fp.startswith('chr')]
    global chrom_dict
    '''
    if bootstrap == True:
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom] = pd.DataFrame()
            ld_proxy_chrom(chrom, query_snps, corr_thresh, window, pop, ld_dir)
    else:
    '''
    manager = mp.Manager()
    chrom_dict = manager.dict()
    for chrom in chrom_list:
        chrom_dict[chrom] = pd.DataFrame()
    with mp.Pool(16) as pool:
        pool.starmap(ld_proxy_chrom,
                     zip(chrom_list,
                        repeat(query_snps),
                        #repeat(chrom_dict),
                        repeat(corr_thresh),
                        repeat(window),
                        repeat(pop),
                        repeat(ld_dir)))
    df = []
    for chrom in chrom_dict.keys():
        df.append(chrom_dict[chrom])
    if len(df) == 0:
        df = pd.DataFrame()
    else:
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


def ld_proxy_chrom_(chrom, query_snps, corr_thresh, window, pop, ld_dir):
    q_snps = query_snps[query_snps['chrom'] == chrom[3:]]['rsid'].tolist()
    if q_snps == []:
        return
    snps = []
    db = create_engine(f'sqlite:///{ld_dir}/{chrom}.db', 
                       echo=False, poolclass=NullPool)
    conn = db.connect()
    metadata = sqlalchemy.MetaData()
    ld = sqlalchemy.Table('ld', metadata, autoload=True, autoload_with=db)
    query = (sqlalchemy
             .select([ld])
             .where(sqlalchemy.and_(ld.columns.rsidq.in_(q_snps),
                                    ld.columns.corr >= corr_thresh))
    )
    res_proxy = conn.execute(query)
    more = True
    while more:
        res = res_proxy.fetchmany(100000)
        if res == []:
            more = False
        else:
            res_df = pd.DataFrame(res, columns=res[0].keys())
            res_df = res_df[(abs(res_df['posq']-res_df['post']) <= window)]
            snps.append(res_df)
    if len(snps) == 0:
        return 
    snps = pd.concat(snps)
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

    
def ld_proxy_(query_snps, corr_thresh, window, pop, ld_dir, logger):
    query_snps =  [snp.strip() for snp in query_snps]
    snps, missed = ns.rsids2pos(query_snps,
                                os.path.join(os.path.dirname(__file__), 'data/snps/'))
    ld_dir = os.path.join(ld_dir, pop)
    chrom_list = [fp.split('.')[0] for fp in os.listdir(ld_dir) if fp.startswith('chr')]
    global chrom_dict
    '''
    if bootstrap == True:
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom] = pd.DataFrame()
            ld_proxy_chrom(chrom, query_snps, corr_thresh, window, pop, ld_dir)
    else:
    '''
    manager = mp.Manager()
    chrom_dict = manager.dict()
    for chrom in chrom_list:
        chrom_dict[chrom] = pd.DataFrame()
    with mp.Pool(16) as pool:
        pool.starmap(ld_proxy_chrom_,
                     zip(chrom_list,
                        repeat(snps),
                        #repeat(chrom_dict),
                        repeat(corr_thresh),
                        repeat(window),
                        repeat(pop),
                        repeat(ld_dir)))
    df = []
    for chrom in chrom_dict.keys():
        df.append(chrom_dict[chrom])
    df = pd.concat(df)
    query_snps_df = (
        pd.DataFrame(
            {'chromq': '',
             'posq': '',
             'rsidq': query_snps,
             'chromt': '',
             'post': '',
             'rsidt': query_snps,
             'corr': 1,
             'dprime': 1
            })
           .merge(snps, how='left', left_on='rsidq', right_on='rsid')
           .assign(#posq = lambda x: x.start,
                   chromq = lambda x: x.chrom,
                   #post = lambda x: x.start,
                   chromt = lambda x: x.chrom)
           .drop(columns=['chrom', 'start', 'rsid'])
    )
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
        help='The r-squared correlation threshold to use. Default = 0.8')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched. Default = 5000')
    parser.add_argument('-p', '--population', default='EUR',
                        choices=['EUR'],
                        help='The ancestral population in which the LD is calculated')
    parser.add_argument(
        '--ld-dir', default=os.path.join(os.path.dirname(__file__), 'data/ld/dbs/super_pop/'),
        help='Directory containing LD database.')
    return parser.parse_args()


if __name__=='__main__':
    start_time = time()
    args = parse_args()
    outdir = os.path.dirname(args.output)
    os.makedirs(outdir, exist_ok=True)
    logger = logger.Logger(logfile=os.path.join(outdir, 'ld_proxy.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    snps = parse_snps(args.snps, logger)
    #df = ld_proxy_(snps['snp'].tolist(), args.correlation_threshold,
    #              args.window, args.population, args.ld_dir, logger)
    df = ld_proxy(snps['snp'].tolist(), args.correlation_threshold,
                   args.window, args.population, args.ld_dir, logger, True)
    write_results(df, args.output)
    logger.write('Total time elasped: {:.2f} mins.'.format(
        (time()-start_time)/60))

