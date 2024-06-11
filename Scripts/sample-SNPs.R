# Script to subset genotypes that are significantly associated with climate as per the GEA
# Needed to speed up linear model genomic prediction models using genotype data

library(data.table)
library(Matrix)

dir = '/group/runciegrp2/Projects/SeeD/'
genotype_dir = paste0(dir, 'Genetic_data/Imputed_V4/genotypes_by_chromosome/')
output_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/sampled_genotypes_by_chromosome/'

data = as.data.frame(fread(file.path(dir, 'Env_data/GEA-climate-invnormtransformed.csv'),data.table = F))
SampleIDs = data[, 1]
non_dup_acc = read.csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/selected_genotypeIDs.csv')
colnames(non_dup_acc) = c('Unique.ID', 'SampleID')

final_genotype = data.frame(matrix(, nrow= 3511, ncol = 0))

for(chr in 1:10){
  mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
  # rownames(mat) = mat[,1]
  mat = as.matrix(mat[,-1])
  mat = mat[match(SampleIDs,non_dup_acc$Unique.ID),]
  rownames(mat) = SampleIDs

  mean_notNA = colMeans(apply(mat,2,function(x) x %in% 0:2))
  maf = .5 - abs(.5-colMeans(mat)/2)
  select_cols = maf > 0.01 & mean_notNA>0.99
  filtered_mat = mat[,select_cols]
  sampled_mat = filtered_mat[,sample(ncol(filtered_mat), size = floor(ncol(filtered_mat)/10) , replace = F)]

  final_genotype = cbind(final_genotype, sampled_mat)
}

row.names(final_genotype) = SampleIDs
write.csv(final_genotype, file.path(output_dir, 'sampled_SNPs_forPCA_noNA.csv'))