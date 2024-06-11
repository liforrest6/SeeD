library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(spatialsample)
library(sf)
library(SpatialKDE)
library(sp)
library(tmap)
library(rsample)
library(foreach)

run = commandArgs(t=T)
rep = as.numeric(run[1])
dir = '/group/runciegrp2/Projects/SeeD/'

print(rep)

K = fread('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

env_data = fread(paste0(dir, 'Env_data/GEA-climate-invnormtransformed.csv'), data.table = F)
geno_info = read.csv(file.path(dir, 'Phenotype_data/selected_genotypeIDs.csv'))
passport_data = readxl::read_excel(file.path(dir, 'Env_data/Passport/Maize climate data nov 2022 AEZ elevation by growing season and flowering period.xlsx'))
dna_to_tc_gid = as.data.frame(readxl::read_xlsx(file.path(dir, 'Phenotype_data/Hearne_data/dataset/Blups_01.xlsx'),sheet = 'DNA to TC GID'))


geno_info$AccID = dna_to_tc_gid[match(geno_info$Sample,dna_to_tc_gid$`Sample ID`),2]
geno_info$LatNew = passport_data$LatNew[match(geno_info$AccID,passport_data$GID)]
geno_info$LongNew = passport_data$LongNew[match(geno_info$AccID,passport_data$GID)]

geno_info = geno_info[!is.na(geno_info$LatNew) & geno_info$V1 %in% colnames(K) & geno_info$V1 %in% env_data$Unique.ID,]
K = K[geno_info$V1,geno_info$V1]
env_data = env_data[match(geno_info$V1,env_data$Unique.ID),]

geno_info = merge(geno_info,env_data,by.x = 'V1',by.y = 'Unique.ID')

geno_sf = st_as_sf(geno_info,coords = c('LongNew','LatNew'))
# need to set coordinate system
st_crs(geno_sf) <- 4326

# clusters = spatial_clustering_cv(geno_sf,v=10,buffer = 4e5)
# clusters_spatial = clusters

clusters = rsample::vfold_cv(geno_sf,v=10)

folds_list = list()
for(i in 1:10){
	z = data.frame(x = assessment(clusters$splits[[i]])$V1, sprintf('Fold%d', i))
	print(nrow(z))
	folds_list[[i]] = z
}



all_folds = do.call('rbind', folds_list)
colnames(all_folds) = c('accession', 'fold')
write.csv(all_folds, file.path(dir, sprintf('Analyses/MegaLMM_output/dataset_clim-spatial_prediction-rep_%02d/testing_folds_assignment.csv', rep)), row.names = F, quote = F)


pdf(file.path(dir, sprintf('Analyses/Spatial/spatial_cluster_map_rep_%02d.pdf', rep)))
autoplot(clusters)
dev.off()


