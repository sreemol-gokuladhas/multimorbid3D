#! /usr/bin/env python
import pandas as pd
import os
import sys
from sqlalchemy import create_engine
from tqdm import tqdm


merged_arch = 'data/snps/human_9606_b151_GRCh38p7_RsMergeArch.pairs.bcp.gz'
out_dir = 'data/snps/dbs/'
os.makedirs(out_dir, exist_ok=True)



db_url = f'sqlite:///{out_dir}/human_9606_b151_GRCh38p7_RsMergeArch.bcp.db'
desc = f'Building merged snps table...'
bar_format = '{desc}: {n_fmt} {unit}'
t = tqdm(total=0, unit='entries', desc=desc, disable=False,
         bar_format=bar_format)
table = 'merged_arch'
chunksize = 200000
db = create_engine(db_url, echo=False)
idx = 1
for df in pd.read_csv(merged_arch, sep='\t', chunksize=chunksize, 
                      header=None, names=['new', 'old']):
    if_exists = 'replace' if t.total == 0 else 'append'
    df.to_sql(table, con=db, if_exists=if_exists, index=False)
    idx += len(df)
    t.total += len(df)
    t.update(len(df))
    t.close()
db.execute('''CREATE INDEX idx_{} ON {} (new)'''.format(
    '{}_{}'.format(table, 'new'), table))
db.execute('''CREATE INDEX idx_{} ON {} (old)'''.format(
    '{}_{}'.format(table, 'old'), table))
    
