#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/LocalAdaptation_run_%A_%a.txt
#SBATCH -e /home/fli21/slurm-log/LocalAdaptation_run_%A_%a.txt
#SBATCH -J localadaptation
#SBATCH -t 24:00:00
#SBATCH --mem 10GB
#SBATCH -n 8
#SBATCH --array=1-4
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu



module load R/4.3.3


# Rscript /group/runciegrp2/Projects/SeeD/Scripts/Optimal_environment_model_unconstrained.R 3 $SLURM_ARRAY_TASK_ID

Rscript /group/runciegrp2/Projects/SeeD/Scripts/Optimal_environment_model_independent.R 3 $SLURM_ARRAY_TASK_ID


