#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/subsetSnps-%j.txt
#SBATCH -e /home/fli21/slurm-log/subsetSnps-%j.txt
#SBATCH -J subset-SNPs
#SBATCH -t 24:00:00
#SBATCH --mem 16GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

echo 'This is subsetting SNPs for genotype analysis.'
module load R/4.2.2

dir='/group/runciegrp2/Projects/SeeD/'

## subset for enriched SNPs
# Rscript ${dir}/Scripts/subset-SNPs.R /Analyses/GEA_output/multivariate_results/clumped/GEA_clumped_SNPs_list.csv GEA_clumped_SNPs_genotype.csv
## subset for matching SNPs
# Rscript ${dir}/Scripts/subset-SNPs.R /Analyses/GEA_output/multivariate_results/clumped/GEA_matching_sampled_SNPs.txt GEA_matching_SNPs_genotype.csv


## subset for enriched SNPs
# Rscript ${dir}/Scripts/subset-SNPs.R /Analyses/GEA_output/multivariate_results_unstructured/clumped/GEA_clumped_SNPs_list.csv GEA_unstructured_clumped_SNPs_genotype.csv
## subset for matching SNPs
Rscript ${dir}/Scripts/subset-SNPs.R /Analyses/GEA_output/multivariate_results_unstructured/clumped/unstructured_gea_matching_sampled_SNPs.txt GEA_unstructured_matching_SNPs_genotype.csv

## generate random sample of high coverage SNPs for PCA analysis
# Rscript ${dir}/Scripts/sample-SNPs.R

