#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/SeeD/Analyses/single_trait_SKAT
#SBATCH -o /group/runciegrp2/Projects/SeeD/Analyses/logs/out-ST_SKAT-%j_%A_%a.out
#SBATCH -e /group/runciegrp2/Projects/SeeD/Analyses/logs/error-ST_SKAT-%j_%A_%a.err
#SBATCH -J ST_SKAT


module load R

Rscript ../../Scripts/single_trait_GWAS_SKAT.R ${SLURM_ARRAY_TASK_ID} > ../../Analyses/logs/single_trait_GWAS_SKAT_${SLURM_ARRAY_TASK_ID}.Rout