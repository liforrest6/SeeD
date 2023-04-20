#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/JointGWAS
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-%j_%A_%a.err
#SBATCH -J run_Cov


module load R

Rscript ../../Scripts/model_covariance_MegaLMM.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/model_covariance_MegaLMM_${SLURM_ARRAY_TASK_ID}.Rout