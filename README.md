# Multimorbid3D

Predict multimorbid conditions of a disease

## Installation

Clone the repository using the following command:

```
https://github.com/Genome3d/multimorbid3D.git
cd multimorbid3D

```

## Setup

1. Create a conda environment from [environment.yaml](https://github.com/Genome3d/multimorbid3D/blob/main/environment.yaml) file

```
conda env create --prefix ./multimorbid3D --file environment.yaml

```

2. Activate the environment 

```
conda activate multimorbid3D/

```

## Usage

```
(/mnt/projects/multimorbid3D/env) sgok603@liggen3dprd03:/mnt/projects/multimorbid3D$ python comorbid.py -h
usage: comorbid.py [-h] [-g GENES [GENES ...]] [-s SNPS [SNPS ...]] [--trait TRAIT] [--pmid PMID] --grn-dir GRN_DIR [--gwas GWAS] -o OUTPUT_DIR [-l LEVELS] [-p {string,proper} [{string,proper} ...]]
                   [--string-score STRING_SCORE] [--bootstrap] [--bootstraps BOOTSTRAPS] [--keep-bootstraps] [--non-spatial] [--non-spatial-dir NON_SPATIAL_DIR] [--snp-ref-dir SNP_REF_DIR]
                   [--gene-ref-dir GENE_REF_DIR] [--ld] [-c CORRELATION_THRESHOLD] [-w WINDOW] [--population {EUR}] [--ld-dir LD_DIR]

Identify multimorbid traits based on eQTL associations and protein-protein interactions.

optional arguments:
  -h, --help            show this help message and exit
  -g GENES [GENES ...], --genes GENES [GENES ...]
                        A space-separated list of gene symbols or filepath to a file containing gene symbols in the 'gene' column.
  -s SNPS [SNPS ...], --snps SNPS [SNPS ...]
                        A space-separated list of SNP rsIDs or filepath to a file containing SNP rsids in the 'snp' column.
  --trait TRAIT         GWAS trait to query. Note: this flag is mutually exclusive with the --snps and --pmid flags
  --pmid PMID           PubMed ID of the GWAS to query. Note: this flag is mutually exclusive with the --snps and --trait flag
  --grn-dir GRN_DIR     Directory containing tissue gene regulatory network. The subdirectories should contain significant_eqtls.txt for each chromosome.
  --gwas GWAS           Filepath to GWAS associations. Default: Associations from the GWAS Catalog (https://www.ebi.ac.uk/gwas/api/search/downloads/full)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to write results.
  -l LEVELS, --levels LEVELS
                        Path length (i.e. number of nodes) to query. Default = 1
  -p {string,proper} [{string,proper} ...], --ppin {string,proper} [{string,proper} ...]
                        The protein-protein-interaction database(s) to use. Default: ['string', 'proper']
  --string-score STRING_SCORE
                        Cut-off score for STRING interactions. Default = 0.7
  --bootstrap           Perform a bootstrap. Default = False
  --bootstraps BOOTSTRAPS
                        Number of bootstrap datasets. Default: 1000
  --keep-bootstraps     Keep bootstrap results. Default: False
  --non-spatial         Include non-spatial eQTLs. Default = False
  --non-spatial-dir NON_SPATIAL_DIR
                        Filepath to non-spatial eQTLs.
  --snp-ref-dir SNP_REF_DIR
                        Filepath to SNP BED databases.
  --gene-ref-dir GENE_REF_DIR
                        Filepath to gene BED.
  --ld                  Include LD SNPs in identifying eQTLs and GWAS traits. Default = False
  -c CORRELATION_THRESHOLD, --correlation-threshold CORRELATION_THRESHOLD
                        The r-squared correlation threshold to use.
  -w WINDOW, --window WINDOW
                        The genomic window (+ or - in bases) within which proxies are searched. Default = 5000
  --population {EUR}    The ancestral population in which the LD is calculated. Default = "EUR"
  --ld-dir LD_DIR       Directory containing LD database.
```
