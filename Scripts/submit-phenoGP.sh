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
module load R/4.2.2

dir='/group/runciegrp2/Projects/SeeD/'

## parameters = transform traitN ntree output mex

# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid
# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid _mex

# Rscript ${dir}/Scripts/subset-enriched-SNPs.R

# for i in {1..5}; do
# 	Rscript ${dir}/Scripts/phenotypic_GP.R 1 ${i} 1000 deregressed_blups_nostress
# done

Rscript ${dir}/Scripts/phenotypic_GP.R 1 $SLURM_ARRAY_TASK_ID 1000 deregressed_blups _nostress