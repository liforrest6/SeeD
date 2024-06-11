# Script to subset genotypes that are significantly associated with climate as per the GEA
# Needed to speed up linear model genomic prediction models using genotype data

library(data.table)

run = commandArgs(t=T)
SNP_file = run[1]
output_file_name = run[2]

dir = '/group/runciegrp2/Projects/SeeD/'
genotype_dir = paste0(dir, 'Genetic_data/Imputed_V4/genotypes_by_chromosome/')

data = as.data.frame(fread(file.path(dir, 'Env_data/GEA-climate-invnormtransformed.csv'),data.table = F))
# SNP_list = as.data.frame(fread(file.path(dir, 'Analyses/GEA_output/multivariate_results/clumped/GEA_clumped_SNPs_list.csv'), header = F))$V1
SNP_list = as.data.frame(fread(file.path(dir, SNP_file), header = F))$V1
SampleIDs = data[, 1]
non_dup_acc = read.csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/selected_genotypeIDs.csv')

final_genotype = data.frame(matrix(, nrow= 3511, ncol = 0))

for(chr in 1:10){
  mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
  rownames(mat) = mat[,1]
  mat = as.matrix(mat[,-1])
  mat = mat[match(SampleIDs,non_dup_acc$V1),colnames(mat) %in% SNP_list]

  final_genotype = cbind(final_genotype, mat)
}

row.names(final_genotype) = SampleIDs
# write.csv(final_genotype, file.path(genotype_dir, 'GEA_clumped_SNPs_genotype.csv'))
write.csv(final_genotype, file.path(genotype_dir, output_file_name))