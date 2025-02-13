library(rgbif)
library(sf)
library(mapview)
library(tmap)
library(dplyr)


cimmyt = read_excel('/Users/liforrest/Documents/Projects/gea-adaptation/data/dataset/cimmyt_data.xlsx', 
                    sheet = 'climate growing and fl season',
                    guess_max = 12000)

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
genotypesWithIDs = read.csv(here(env_data_dir, 'GenotypesWithIDs.csv'))[c(1,2,6,11)]

acc_coords = merge(merge(finalMat, genotypesWithIDs, by = 'Unique.ID'), cimmyt[c(1,4,5,7,8)], by = 'accid')

acc_points = acc_coords %>% st_as_sf(coords = c('LongNew', 'LatNew'),
                                     crs = '+proj=longlat +datum=WGS84 +ellps=WGS84')
mapview(acc_points)
mex_points = acc_points %>% filter(oiso == 'MEX') %>% filter(LongNew > -106.8 & LongNew < -95.5)
mapview(mex_points)

filt_bbox = sf::st_bbox(c(xmin = -106.8,
                          ymin = 15.5,
                          xmax = -95.5,
                          ymax = 22.8),
                        crs = st_crs('+proj=longlat +datum=WGS84 +ellps=WGS84')) %>% 
  sf::st_as_sfc(.)
find_acc = sf::st_within(mex_points, filt_bbox)
filt_acc = mex_points[which(lengths(find_acc) != 0), ]
mapview(filt_acc)
filt_acc$Unique.ID

filt_df = as.data.frame(filt_acc$Unique.ID)
colnames(filt_df) = c('Unique.ID.Mex')
# write.csv(filt_df, 'Env_data/central_mexican_accessions.csv')

all_coords = cimmyt[1:9] %>% st_as_sf(coords = c('LongNew', 'LatNew'),
                                      crs = '+proj=longlat +datum=WGS84 +ellps=WGS84') 
all_coords = cimmyt[1:9]
terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = -50, right = -33, top = 33), 
                                  zoom = 4, maptype = "terrain-background")

ggmap(terrain_kernels)+
  geom_point(data = all_coords, 
             aes(x = LongNew, y = LatNew, color = AltM),
             pch=16,alpha=0.2,size=1) + 
  scale_color_continuous(name = "Elevation (m)", low = 'black', high = 'red') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 16)) +
  # ggtitle("Geo-coordinates for CIMMYT data", ) +
  xlab("") +
  ylab("")



terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = -40, right = -35, top = 33), 
                                  zoom = 4, maptype = "toner-lite", color = 'bw', messaging = F)




