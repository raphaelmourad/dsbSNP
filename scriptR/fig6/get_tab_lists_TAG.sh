#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --job-name="TAG_test_pval"
#SBATCH --mem=175GB
#SBATCH --output=TAG_test_pval.out

cd /mnt/DATA/MCD/eq_legube/SEB/GWAS/scriptR_npval



#args SNP type lists and number of pval under 10-6 in GWAS
Rscript --vanilla get_tab_lists_TAG.R 'expe.txt' 500
