####################################################################################
# Process climate data for GEA analyses
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################



library(ggplot2)
library(GGally)
library(ggfortify)
library(ggmap)


## read CIMMYT 2023 data
cimmyt = read_excel('/Users/liforrest/Documents/Projects/gea-adaptation/data/dataset/cimmyt_data.xlsx', 
                    sheet = 'climate growing and fl season',
                    guess_max = 12000)
## read elevation data to append
cimmyt_elevation = read_excel('/Users/liforrest/Documents/Projects/gea-adaptation/data/dataset/cimmyt_data.xlsx', sheet = 'elevation')
genotypesWithIDs = read.csv(here(env_data_dir, 'GenotypesWithIDs.csv'))[c(1,2,6,11)]


getTossOutAccessions = function(df) {
  acc_9999 = filter_at(df, vars(starts_with('prec')), any_vars(. == -9999)) %>% pull(accid)
  no_elevation = c('CIMMYTMA-002934', 'CIMMYTMA-001393', 'CIMMYTMA-000254')
  result = append(acc_9999, no_elevation)
  return(result)
}

tossOutAccessions = function(df, list) {
  #drop -9999 values like isla mujer
  result = df %>% filter(!accid %in% list)
  return(result)
}

## toss out 5 accessions, total = 3511 accessions
tossOutList = getTossOutAccessions(cimmyt)

cimmyt_grow = merge(genotypesWithIDs, 
                    cimmyt, by.x = 'accid',
                    by.y = 'accid') %>% 
  filter(!accid %in% tossOutList) %>% 
  rename('prec_01'='prec_01...10',
         'prec_02'='prec_02...11',
         'prec_03'='prec_03...12',
         'prec_04'='prec_04...13',
         'prec_05'='prec_05...14', 
         'prec_06'='prec_06...15',
         'prec_07'='prec_07...16',
         'prec_08'='prec_08...17',
         'prec_09'='prec_09...18',
         'prec_10'='prec_10...19',
         'prec_11'='prec_11...20',
         'prec_12'='prec_12...21',
         'gs_prec_01'='prec_01...39',
         'gs_prec_02'='prec_02...40',
         'gs_prec_03'='prec_03...41',
         'gs_prec_04'='prec_04...42',
         'gs_prec_05'='prec_05...43',
         'gs_prec_06'='prec_06...44',
         'gs_prec_07'='prec_07...45',
         'gs_prec_08'='prec_08...46',
         'gs_prec_09'='prec_09...47',
         'gs_prec_10'='prec_10...48',
         'gs_prec_11'='prec_11...49',
         'gs_prec_12'='prec_12...50')

## calculate aggregate data for growing season months (default except for precipitation)
cimmyt_grow$tmin = cimmyt_grow %>% dplyr::select(starts_with('tmin')) %>% apply(., 1, FUN = min, na.rm = T)
cimmyt_grow$tmax = cimmyt_grow %>% dplyr::select(starts_with('tmax')) %>% apply(., 1, FUN = max, na.rm = T)
cimmyt_grow$trange = cimmyt_grow$tmax - cimmyt_grow$tmin
cimmyt_grow$precipTot = cimmyt_grow %>% dplyr::select(starts_with('gs_prec')) %>% apply(., 1, FUN = sum, na.rm = T)
cimmyt_grow$aridityMean = cimmyt_grow %>% dplyr::select(starts_with('aridity')) %>% apply(., 1, FUN = mean, na.rm = T)
cimmyt_grow$rhMean = cimmyt_grow %>% dplyr::select(starts_with('rhumav')) %>% apply(., 1, FUN = mean, na.rm = T)

cimmyt_elevation = tossOutAccessions(cimmyt_elevation, tossOutList)
cimmyt_grow$elevation = cimmyt_elevation[match(cimmyt_grow$accid, cimmyt_elevation$accid),]$elevation

## dataframe using summary statistics that were already in CIMMYT dataset
cimmyt_established = cimmyt_grow[,c('accid', 'GID', 'LatNew', 'LongNew',
                                    'Cropping season total pptn = 6m max',
                                    'Sum pptn over growing season',
                                    'Pooled Mean pptn growing season',
                                    'Sum pptn over fl',
                                    'Pooled Mean pptn fl',
                                    'Pooled Mean Tmax',
                                    'Pooled Mean fl Tmax',
                                    'Pooled Mean Tmin',
                                    'Pooled Mean pptn Tmin',
                                    'Pooled monthly mean aridity index growing season (scaled)',
                                    'Sum aridity over fl',
                                    'Pooled Mean aridity index fl (scaled)',
                                    'Pooled mean rhum growing season',
                                    'Pooled Mean rhum Tmin')]

## check correlation between data
initial_columns = c('Sample.ID', 'LongNew', 'LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
final_columns = c('Unique.ID', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')

invNormTransform = function(x) {
  qnorm((rank(x,na.last="keep")-3/8)/sum(!is.na(x)))
}

finalMat = cimmyt_grow[final_columns]
finalMat_transform = finalMat
finalMat_transform[2:8] = lapply(finalMat[2:8], FUN = invNormTransform)

if(!file.exists(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))) {
  write.csv(finalMat_transform, here(env_data_dir, 'GEA-climate-invnormtransformed.csv'), row.names = F)
}
