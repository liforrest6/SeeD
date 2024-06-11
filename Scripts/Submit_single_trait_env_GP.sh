#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/GenomicPrediction
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-ST_GP-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-ST_GP-%j_%A_%a.err
#SBATCH -J ST_GP


module load R

Rscript ../../Scripts/single_trait_env_GP.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/single_trait_env_GP_${SLURM_ARRAY_TASK_ID}.Rout