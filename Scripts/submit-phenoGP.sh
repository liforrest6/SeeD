#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/phenoGP-%j.txt
#SBATCH -e /home/fli21/slurm-log/phenoGP-%j.txt
#SBATCH -J phenoGP
#SBATCH -t 24:00:00
#SBATCH --mem 8GB
#SBATCH -n 8
#SBATCH --array=1-5
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running a phenotypic GP'
module load R/4.3.3

dir='/group/runciegrp2/Projects/SeeD/'

## parameters = transform traitN ntree output nostress

# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid
# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid mex

# Rscript ${dir}/Scripts/subset-enriched-SNPs.R

## iterative before using slurm array
# for i in {1..5}; do
# 	Rscript ${dir}/Scripts/phenotypic_GP.R 1 ${i} 1000 deregressed_blups_nostress
# done

## included _filterTrial trait to take out tester:trial combinations with too low sample size
# Rscript ${dir}/Scripts/phenotypic_GP.R 1 $SLURM_ARRAY_TASK_ID 1000 deregressed_blups nostress 

## testing use of unstructured envGWAS SNPs to see how it compares to structured GEA
Rscript ${dir}/Scripts/phenotypic_GP_unstructured.R 1 $SLURM_ARRAY_TASK_ID 1000 deregressed_blups nostress 