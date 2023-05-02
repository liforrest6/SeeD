##########################################
# Process genotype set for GEA analyses
#
# Author: Forrest Li
# Script for fill season
##########################################

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

