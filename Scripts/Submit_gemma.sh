#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/JointGWAS
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-gemma-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-gemma-%j_%A_%a.err
#SBATCH -J gemma


module load R

Rscript ../../Scripts/single_trait_GWAS_GEMMA.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/single_trait_GWAS_GEMMA_${SLURM_ARRAY_TASK_ID}.Rout