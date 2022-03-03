#! /usr/bin/env python
import pandas as pd
import os
import sys
from sqlalchemy import create_engine
from tqdm import tqdm


ld_dir = 'data/ld/super_pop/EUR'
out_dir = 'data/ld/dbs/super_pop/EUR'
os.makedirs(out_dir, exist_ok=True)


for chrom_fp in os.listdir(ld_dir):
    chrom = os.path.basename(chrom_fp)
    if chrom == 'chr16':
        continue
    db_url = f'sqlite:///{out_dir}/{chrom}.db'
    desc = f'Building {chrom} LD table...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='assoc', desc=desc, disable=False,
             bar_format=bar_format)
    table = 'ld'
    chunksize = 200000
    cols = ['chromq', 'posq', 'rsidq', 'chromt', 'post', 'rsidt', 'corr', 'dprime']
    db = create_engine(db_url, echo=False)
    idx = 1
    chrom_fp = os.path.join(ld_dir, chrom_fp)
    for df in pd.read_csv(chrom_fp, sep=',', chunksize=chunksize, skiprows=1,
                          header=None, names=list(range(22)), usecols=[0,1,4,8,9,12,20,21]):
        df.columns = cols
        if_exists = 'replace' if t.total == 0 else 'append'
        df.to_sql(table, con=db, if_exists=if_exists, index=False)
        idx += len(df)
        t.total += len(df)
        t.update(len(df))
    t.close()
    db.execute('''CREATE INDEX idx_{} ON {} (rsidq)'''.format(
        '{}_{}'.format(table, 'rsidq'), table))
    
