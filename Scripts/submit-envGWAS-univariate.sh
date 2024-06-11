#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/envGWAS-univariate-%j.txt
#SBATCH -e /home/fli21/slurm-log/envGWAS-univariate-%j.txt
#SBATCH -J envGWAS-uni
#SBATCH -t 4:00:00
#SBATCH --mem 40GB
#SBATCH -n 2
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running a univariate envGWAS.'
module load R/4.2.2
module load gemma/0.98.5

dir='/group/runciegrp2/Projects/SeeD/'

Rscript ${dir}/Scripts/envGWAS-univariate.R 4
# > ${dir}Analyses/logs/envGWAS-univariate_${SLURM_ARRAY_TASK_ID}.Rout