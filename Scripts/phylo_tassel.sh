#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/phylo-tassel-%j.txt
#SBATCH -e /home/fli21/slurm-log/phylo-tassel-%j.txt
#SBATCH -J phylo-tassel
#SBATCH -t 24:00:00
#SBATCH --mem 128GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

module load bcftools
module load jdk
module load tassel

seed_directory='/group/runciegrp2/Projects/SeeD/'


# echo Starting zipping job...
# bgzip -c ${seed_directory}/Genetic_data/Imputed_V4/fullMaizeGBS.vcf > ${seed_directory}/Genetic_data/Imputed_V4/fullMaizeGBS.vcf.gz
# echo Done zipping and starting indexing...
# tabix -p vcf ${seed_directory}/Genetic_data/Imputed_V4/fullMaizeGBS.vcf.gz
echo Done indexing and starting filtering...
bcftools view --force-samples --samples-file ${seed_directory}/Env_data/finalGEA_accession_list.txt ${seed_directory}/Genetic_data/Imputed_V4/fullMaizeGBS.vcf.gz > ${seed_directory}/Genetic_data/Imputed_V4/maizeGBS_GEA.subset.vcf
echo Finished filtering

echo Starting PCA...
run_pipeline.pl -Xmx128g -fork1 -importGuess ${seed_directory}/Genetic_data/Imputed_V4/maizeGBS_GEA.subset.vcf \
-PrincipalComponentsPlugin -covariance true -endPlugin -export ${seed_directory}/Analyses/tassel_output/maizeGBS_GEA.subset.pca -runfork1

echo Starting tree creation...
run_pipeline.pl -Xmx128g -fork1 -importGuess ${seed_directory}/Genetic_data/Imputed_V4/maizeGBS_GEA.subset.vcf -tree Neighbor -export ${seed_directory}/Analyses/tassel_output/maizeGBS_GEA.subset.tree -runfork1 
