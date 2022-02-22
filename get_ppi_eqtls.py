#!/usr/bin/env python
import pandas as pd
import sys
import os
import argparse

def parse_ppin(ppin_fp):
    #print('Parsing input...')
    if os.path.isfile(ppin_fp): 
        df = pd.read_csv(ppin_fp, sep='\t')
        return df

def get_snps(ppin, eqtl_dir, output_dir):
    #print('Identifying node eQTLs...')
    levels = [col for col in ppin.columns if col.startswith('level')]
    res = {}
    for level in levels:
        res[level] = []
    for chrom in os.listdir(eqtl_dir):
        if not chrom.startswith('chr'):
            continue
        #print(f'Querying {chrom}...')
        chrom_res = []
        eqtls = pd.read_csv(
            os.path.join(eqtl_dir, chrom, 'significant_eqtls.txt'),
            sep='\t', usecols=['snp', 'gene'])
        for level in levels:
            res[level].append(eqtls[eqtls['gene'].isin(ppin[level])].drop_duplicates())
    write_results(res, output_dir)
    return res


def write_results(res, out):
    #print('Writing results...')
    os.makedirs(out, exist_ok=True)
    for level in res:
        level_fp = os.path.join(out, f'{level}.txt')
        if len(res[level]) == 0:
            continue
        df = pd.concat(res[level])
        df.to_csv(level_fp, sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Find eQTLs of proteins.')
    parser.add_argument(
        '-p', '--ppin', required=True,
        help='''Filepath to PPIN file of levels of interactions.''' )
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
         '--grn-dir', required=True, help='Filepath to tissue eQTL map.')

    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    ppin = parse_ppin(args.ppin)
    res = get_snps(ppin, args.grn_dir, args.output_dir)

