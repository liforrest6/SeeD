#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/randomForest-%j.txt
#SBATCH -e /home/fli21/slurm-log/randomForest-%j.txt
#SBATCH -J RF
#SBATCH -t 24:00:00
#SBATCH --mem 8GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running a random forest.'
module load R/4.2.2

dir='/group/runciegrp2/Projects/SeeD/'


# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid
# Rscript ${dir}/Scripts/randomforest_env_GP.R 1 1 1000 resid _mex

# Rscript ${dir}/Scripts/subset-enriched-SNPs.R

Rscript ${dir}/Scripts/phenotypic_GP.R 1 1 1000 nonG