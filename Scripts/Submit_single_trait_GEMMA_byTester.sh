#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/JointGWAS
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-single_trait_GEMMA_byTester-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-single_trait_GEMMA_byTester-%j_%A_%a.err
#SBATCH -J ST_G_byT


module load R

Rscript ../../Scripts/single_trait_GWAS_GEMMA_byTester.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/single_trait_GWAS_GEMMA_byTester_${SLURM_ARRAY_TASK_ID}.Rout