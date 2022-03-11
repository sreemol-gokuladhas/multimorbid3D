#!/usr/bin/env python
import pandas as pd
import sys
import os
import argparse

import query_grn

def parse_ppin(ppin_dir):
    #print('Parsing input...')
    if os.path.isdir(ppin_dir):
        gene_fps = sorted([level for level in os.listdir(ppin_dir)
                     if level.endswith('_genes.txt')])
        gene_list = []
        for level in gene_fps:
            df = pd.read_csv(os.path.join(ppin_dir, level), sep='\t')
            gene_list.append(df['gene'].tolist())
        return gene_list
    else:
        sys.exit('ppin-dir entered not a directory.')

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
        '--ppin-dir', required=True,
        help='''Directory containing files of PPIN level genes.''' )
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
         '--grn-dir', required=True, help='Filepath to tissue eQTL map.')
    parser.add_argument(
        '--non-spatial', action='store_true', default=False,
        help='Include non-spatial eQTLs.')
    parser.add_argument(
        '--non-spatial-dir', default='data/GTEx/', help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '--snp-ref-dir', default='data/snps/', help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default='data/genes/', help='Filepath to gene BED.')
    parser.add_argument(
        '--ld', action='store_true', default=False,
        help='Include LD SNPs in identifying eQTLs and GWAS traits.')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched.')
    parser.add_argument('-p', '--population', default='EUR',
                        choices=['EUR'],
                        help='The ancestral population in which the LD is calculated.')
    parser.add_argument('--ld-dir', default='data/ld/dbs/super_pop/',
                        help='Directory containing LD database.') 
    
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    gene_list = parse_ppin(args.ppin_dir)
    grn = query_grn.parse_grn(args.grn_dir)
    query_grn.get_gene_eqtls(
        gene_list, grn, args.output_dir,
        args.non_spatial, args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir)
