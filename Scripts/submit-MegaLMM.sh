#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/MegaLMM-prediction-%j.txt
#SBATCH -e /home/fli21/slurm-log/MegaLMM-prediction-%j.txt
#SBATCH -J MegaLMM-prediction
#SBATCH -t 24:00:00
#SBATCH --mem 128GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running MegaLMM for new CIMMYT data for 7 variables.'
module load R/4.2.2

dir='/group/runciegrp2/Projects/SeeD/'

## generates MegaLMM
# Rscript ${dir}/Scripts/envGWAS-MegaLMM.R 2 clim
## generates prediction for MegaLMM with one withheld set
# Rscript ${dir}/Scripts/envGWAS-MegaLMM-prediction.R 1 clim
## performs k-fold prediction
# Rscript ${dir}/Scripts/envGWAS-MegaLMM-prediction-kfold.R 5 clim



## performs spatial
Rscript ${dir}/Scripts/envGWAS-MegaLMM-prediction-spatial.R 2 clim

