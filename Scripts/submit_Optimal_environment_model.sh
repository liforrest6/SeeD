#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/LocalAdaptation
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/LocalAdaptation/logs/run_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/LocalAdaptation/logs/run_%A_%a.err
#SBATCH -J Opt



module load R


Rscript /group/runciegrp2/Projects/SeeD/Scripts/Optimal_environment_model.R 1 1 $SLURM_ARRAY_TASK_ID > logs/Optimal_environment_model_1_1_$SLURM_ARRAY_TASK_ID.Rout

