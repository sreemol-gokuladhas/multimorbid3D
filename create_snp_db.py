#! /usr/bin/env python
import pandas as pd
import os
import sys
from sqlalchemy import create_engine
from tqdm import tqdm


snp_dir = 'data/snps/human_9606_b151_GRCh38p7'
out_dir = 'data/human_9606_b151_GRCh38p7'
os.makedirs(out_dir, exist_ok=True)


for chrom_fp in os.listdir(snp_dir):
    chrom = chrom_fp.split('.')[0]
    db_url = f'sqlite:///{out_dir}/{chrom}.db'
    desc = f'Building {chrom} SNPs table...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='variants', desc=desc, disable=False,
             bar_format=bar_format)
    table = 'snps'
    chunksize = 200000
    cols = ['chr', 'frag_id', 'id', 'variant_id']
    db = create_engine(db_url, echo=False)
    idx = 1
    chrom_fp = os.path.join(snp_dir, chrom_fp)
    for df in pd.read_csv(chrom_fp, sep='\t', chunksize=chunksize, skiprows=1,
                          header=None):
        df.columns = ['chrom', 'start', 'end', 'rsid', '1', '2']
        use_cols = ['chrom', 'start', 'rsid']
        df['chrom'] = df['chrom'].str[3:]
        if_exists = 'replace' if t.total == 0 else 'append'
        df[use_cols].to_sql(table, con=db, if_exists=if_exists, index=False)
        idx += len(df)
        t.total += len(df)
        t.update(len(df))
    t.close()
    db.execute('''CREATE INDEX idx_{} ON {} (rsid)'''.format(
        '{}_{}'.format(table, 'rsid'), table))
    db.execute('''CREATE INDEX idx_{} ON {} (chrom, start)'''.format(
        '{}_{}'.format(table, 'chrom_start'), table))
    
