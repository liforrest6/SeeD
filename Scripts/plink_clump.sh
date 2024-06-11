#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/plink-%j.txt
#SBATCH -e /home/fli21/slurm-log/plink-%j.txt
#SBATCH -J plink-clump
#SBATCH -t 24:00:00
#SBATCH --mem 64GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=med2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is running plink LD analysis'
module load plink/1.9-beta7.1

plink --vcf /group/runciegrp2/Projects/SeeD/Genetic_data/Unimputed_V4/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.vcf \
--keep /group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/selected_genotypeIDs.txt \
--r2 --ld-window 1000 --ld-window-kb 100000 -ld-window-r2 0.1 \
--ld-snp-list /group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/GEA_SNP_1e2_forplink.txt \
--out /group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/plink/LD-analysis

