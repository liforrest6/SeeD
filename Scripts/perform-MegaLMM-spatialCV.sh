#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/MegaLMM-prediction-%j.txt
#SBATCH -e /home/fli21/slurm-log/MegaLMM-prediction-%j.txt
#SBATCH -J MegaLMM-spatial
#SBATCH -t 24:00:00
#SBATCH --mem 200GB
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --array=1-10


echo 'This is running MegaLMM cross-validation for spatial folds'
module load R/4.2.2

repnum=3
dir='/group/runciegrp2/Projects/SeeD/'

mkdir -p ${dir}/Analyses/MegaLMM_output/dataset_clim-spatial_prediction-rep_0${repnum}

## use spatial sample to create fold
# Rscript create_spatial_folds.R ${repnum}

## run actual prediction with MegaLMM, run number then fold number
# for fold in 1 .. 10 do {
# 	Rscript spatial-prediction-byFold.R 2 ${fold}
# } &done

## run actual prediction with MegaLMM, run number then fold number
# Rscript spatial-prediction-byFold.R 2 1

Rscript spatial-prediction-byFold.R ${repnum} ${SLURM_ARRAY_TASK_ID}
