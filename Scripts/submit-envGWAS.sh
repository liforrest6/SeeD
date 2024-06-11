#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/envGWAS-multivariate-%j.txt
#SBATCH -e /home/fli21/slurm-log/envGWAS-multivariate-%j.txt
#SBATCH -J envGWAS
#SBATCH -t 24:00:00
#SBATCH --mem 256GB
#SBATCH -n 64
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running a multivariate envGWAS.'
module load R/4.2.2

dir='/group/runciegrp2/Projects/SeeD/'


Rscript ${dir}/Scripts/envGWAS-multivariate.R clim 9 10
# Rscript ${dir}/Scripts/envGWAS-univariate.R clim
# Rscript ${dir}/Scripts/envGWAS-multivariate.R clim_perm