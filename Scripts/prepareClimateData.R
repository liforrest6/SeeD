library(readxl)
library(dplyr)
library(stringr)
library(tidyverse)
library(car)

## read CIMMYT 2023 data
cimmyt = read_excel('/Users/liforrest/Documents/Projects/gea-adaptation/data/dataset/cimmyt_data.xlsx', 
                    sheet = 'climate growing and fl season',
                    guess_max = 12000)
## read elevation data to append
cimmyt_elevation = read_excel('/Users/liforrest/Documents/Projects/gea-adaptation/data/dataset/cimmyt_data.xlsx', sheet = 'elevation')


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

tossOutList = getTossOutAccessions(cimmyt)

cimmyt_grow = tossOutAccessions(cimmyt, tossOutList) %>% 
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



## set values for grain fill season
cimmyt_fill = tossOutAccessions(cimmyt, tossOutList) %>% 
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
         'fill_prec_01'='prec_01...39',
         'fill_prec_02'='prec_02...40',
         'fill_prec_03'='prec_03...41',
         'fill_prec_04'='prec_04...42',
         'fill_prec_05'='prec_05...43',
         'fill_prec_06'='prec_06...44',
         'fill_prec_07'='prec_07...45',
         'fill_prec_08'='prec_08...46',
         'fill_prec_09'='prec_09...47',
         'fill_prec_10'='prec_10...48',
         'fill_prec_11'='prec_11...49',
         'fill_prec_12'='prec_12...50')
## make columns for first, mid, and last month of growing season
cimmyt_fill[c('Fill_first_month', 'Fill_last_month')] = str_split_fixed(cimmyt_fill$`Flowering and grain fill code`, ' ', 2)
find_mid_month = function(x) {
  result = switch(
    x,
    '01' = '02',
    '02' = '03',
    '03' = '04',
    '04' = '05',
    '05' = '06',
    '06' = '07',
    '07' = '08',
    '08' = '09',
    '09' = '10',
    '10' = '11',
    '11' = '12',
    '12' = '01',
  )
  return(result)
}
cimmyt_fill$Fill_mid_month = sapply(cimmyt_fill$Fill_first_month, FUN = find_mid_month)

## for each month, look for rows where grain fill months is false and assign those columns to NA
for(month in 1:12) {
  for(column in c(sprintf('fill_prec_%02d', month),
                  sprintf('tmin_%02d', month),
                  sprintf('tmax_%02d', month),
                  sprintf('aridity%02d', month),
                  sprintf('rhumav%02d', month))) {
    cimmyt_fill[! (cimmyt_fill$Fill_first_month == sprintf('%02d', month) |
                     cimmyt_fill$Fill_mid_month == sprintf('%02d', month) |
                     cimmyt_fill$Fill_last_month == sprintf('%02d', month) ), column] = NA
  }
}

cimmyt_fill$tmin = cimmyt_fill %>% dplyr::select(starts_with('tmin')) %>% apply(., 1, FUN = min, na.rm = T)
cimmyt_fill$tmax = cimmyt_fill %>% dplyr::select(starts_with('tmax')) %>% apply(., 1, FUN = max, na.rm = T)
cimmyt_fill$trange = cimmyt_fill$tmax - cimmyt_fill$tmin
cimmyt_fill$precipTot = cimmyt_fill %>% dplyr::select(starts_with('gs_prec')) %>% apply(., 1, FUN = sum, na.rm = T)
cimmyt_fill$aridityMean = cimmyt_fill %>% dplyr::select(starts_with('aridity')) %>% apply(., 1, FUN = mean, na.rm = T)
cimmyt_fill$rhMean = cimmyt_fill %>% dplyr::select(starts_with('rhumav')) %>% apply(., 1, FUN = mean, na.rm = T)
cimmyt_fill$elevation = tossOutAccessions(cimmyt_elevation, tossOutList)$elevation


library(ggplot2)

ggplot(cimmyt_grow, aes(x = tmin)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = tmax)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = trange)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = precipTot)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = aridityMean)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = rhMean)) +
  geom_histogram()

ggplot(cimmyt_grow, aes(x = elevation)) +
  geom_histogram()


### sanity check map for elevation and locations
library(ggmap)

terrain_kernels <- get_stamenmap( bbox = c(left = -160, bottom = -50, right = -30, top = 33), 
                                  zoom = 4, maptype = "terrain-background")

ggmap(terrain_kernels)+
  geom_point(data = cimmyt_grow, 
             aes(x = LongNew, y = LatNew, color = elevation),
             pch=16,alpha=0.1,size=1) + 
  scale_color_continuous(name = "Elevation (m)", low = 'black', high = 'red') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 16)) +
  ggtitle("Geo-coordinates for CIMMYT data", ) +
  xlab("") +
  ylab("")


## check correlation between data
initial_columns = c('accid', 'GID', 'LongNew', 'LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
library(GGally)
library(ggfortify)

ggplot(cimmyt_grow)

invNormTransform = function(x) {
  qnorm((rank(x,na.last="keep")-3/8)/sum(!is.na(x)))
}

finalMat = cimmyt_grow[initial_columns]
finalMat_transform = finalMat
finalMat_transform[5:11] = lapply(finalMat[5:11], FUN = invNormTransform)

pca_grow = prcomp(finalMat[5:11], scale = T)
autoplot(pca_grow, colour = 'elevation')

pca_grow_transform = prcomp(finalMat_transform[5:11])
autoplot(pca_grow_transform, colour = 'elevation')

ggpairs(finalMat_transform[4:11])
pcor(finalMat_transform[4:11])

ggmap(terrain_kernels)+
  geom_point(data = finalMat, 
             aes(x = LongNew, y = LatNew, color = rhMean),
             pch=16,alpha=0.2,size=1) + 
  scale_color_continuous(name = "Elevation (m)", low = 'black', high = 'red') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 16)) +
  ggtitle("Geo-coordinates for CIMMYT data", ) +
  xlab("") +
  ylab("")


plot(cimmyt_grow$elevation, finalMat$elevation)
plot(cimmyt_grow$tmin, finalMat$tmin)
plot(cimmyt_grow$tmax, finalMat$tmax)
plot(cimmyt_grow$trange, finalMat$trange)
plot(cimmyt_grow$precipTot, finalMat$precipTot)
plot(cimmyt_grow$aridityMean, finalMat$aridityMean)
plot(cimmyt_grow$rhMean, finalMat$rhMean)

library(gridExtra)
grid.arrange()
