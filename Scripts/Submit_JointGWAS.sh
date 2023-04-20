#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/JointGWAS
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-JointGWAS-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-JointGWAS-%j_%A_%a.err
#SBATCH -J run_JointGWAS


module load R

Rscript ../../Scripts/JointGWAS.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/JointGWAS_${SLURM_ARRAY_TASK_ID}.Rout

