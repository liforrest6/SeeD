library(readxl)
library(tidyverse)
library(data.table)
library(dplyr)
library(raster)
library(sp)
library(MASS)

# using TC_info for coordinates (BLUPS data)
# TC_info = as.data.frame(read_xlsx('/group/runciegrp2/Projects/SeeD_GWAS_Gates/Downloaded_data/Hearne_data/dataset/Blups_01.xlsx',sheet = 'DNA to TC GID', col_types = 'text'))
# using Trial_info for coordinates
# TC_info = as.data.frame(read.csv('/home/fli21/dangates-egwas/prepped_data/SEEDGWAS_passport_environ_TC_info.csv'))

TC_info = fread('/group/runciegrp2/Projects/SeeD/Env_data/GEA-climate-nontransformed.csv')
TC_info = rename(TC_info, longitude = LongNew, latitude = LatNew) 
TC_info_points = subset(TC_info,!is.na(longitude+latitude))


# subset_info = TC_info[,1:3] %>% drop_na()
# lats = na.omit(subset_info$latitude)
# lats_fraction = fractions(lats-floor(lats))
# lats_fraction = data.frame(row = which(!is.na(subset_info$latitude)),do.call(rbind,lapply(as.character(lats_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
# subset_info$denominators = lats_fraction$X2
# subset_info[subset_info==0]=1



raster_file = list.files(path = '/group/runciegrp/SharedResources/WorldClim/v2/',pattern = '.tif',full.names = T)

raster_file = raster_file[-c(20:32)]
worldclim_rasters = lapply(raster_file,function(x) raster(x))
rasters_annotation = fread('/group/runciegrp/SharedResources/WorldClim/v2/Rasters_annotation.csv',data.table=F)
names(worldclim_rasters) = rasters_annotation$Name[match(basename(raster_file),rasters_annotation$raster)]
worldclim = stack(worldclim_rasters)

TC_worldclim = raster::extract(worldclim,SpatialPoints(cbind(TC_info_points$longitude,TC_info_points$latitude)))
TC_worldclim = data.frame(TC_info_points[,c('Sample.ID','Unique.ID','longitude','latitude')],TC_worldclim)

write.csv(TC_worldclim, file = '/group/runciegrp2/Projects/SeeD/Env_data/SEEDGWAS_worldclim.csv',row.names = F)
print('done')

## this code looks at variance obtained from geographic variance obtained from looking at uncertainty around geo-coordinates 
# full_adjacents = data.frame()

# for(i in c(-1, 0, 1)){
#   for(j in c(-1, 0, 1)){
#     print(c(i, j))
#     adjacents = subset_info
#     adjacents$latitude = adjacents$latitude + (i/adjacents$denominators)
#     adjacents$longitude = adjacents$longitude + (j/adjacents$denominators)
#     full_adjacents = rbind(full_adjacents, adjacents)
#   }
# }

# also looking at adjacents
# adjacents_worldclim = raster::extract(worldclim,SpatialPoints(cbind(full_adjacents$longitude,full_adjacents$latitude)))
# adjacents_worldclim = data.frame(full_adjacents[,c('SampleID','longitude','latitude')],adjacents_worldclim)


# write.csv(adjacents_worldclim, file = '/home/fli21/gea-adaptation/data/SEEDGWAS_worldclim_variance.csv',row.names = F)