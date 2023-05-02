####################################################################################
# Generate plots
#
# Author: Forrest Li
# Plots from analyses
####################################################################################

source(here::here('config.R'))

manhattan(clim_results, 
          snp = 'SNP',
          bp = 'BP',
          main = 'Multivariate GEA across 7 climate variables, maf > 0.01',
          highlight = top_hits_clumped)

manhattan(clim_results %>% filter(maf > 0.05), 
          snp = 'SNP',
          bp = 'BP',
          main = 'Multivariate GEA across 7 climate variables, maf > 0.05')

manhattan(clim_results %>% filter(CHR == 9 & maf > 0.01), 
          snp = 'SNP',
          bp = 'BP',
          main = 'HSF9, maf > 0.05',
          highlight = clim_results %>% filter(CHR == 9 & BP > 148309510 & BP < 148717470 ) %>% pull(SNP),)

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


terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = 10, right = -80, top = 33), 
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

pca_grow = prcomp(finalMat[5:11], scale = T)
autoplot(pca_grow, colour = 'elevation')

pca_grow_transform = prcomp(finalMat_transform[5:11])
autoplot(pca_grow_transform, colour = 'elevation')

ggpairs(finalMat[4:11])
ggpairs(finalMat_transform[4:11])
pcor(finalMat_transform[4:11])



plot(finalMat_transform$elevation, finalMat$elevation)
plot(finalMat_transform$tmin, finalMat$tmin)
plot(finalMat_transform$tmax, finalMat$tmax)
plot(finalMat_transform$trange, finalMat$trange)
plot(finalMat_transform$precipTot, finalMat$precipTot)
plot(finalMat_transform$aridityMean, finalMat$aridityMean)
plot(finalMat_transform$rhMean, finalMat$rhMean)

ggmap(terrain_kernels)+
  geom_point(data = test %>% filter(oiso %in% c('MEX', 'GUA')), 
             aes(x = LongNew, y = LatNew),
             pch=16,alpha=0.1,size=2, color = 'red') + 
  # scale_color_continuous(name = "Elevation (m)", low = 'black', high = 'red') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 16)) +
  ggtitle("Geo-coordinates for CIMMYT data", ) +
  xlab("") +
  ylab("")





####################################################################################
## plot GEMMA univariate manhattan plots
####################################################################################
gemma_results = list(tmin_results, tmax_results, trange_results, precipTot_results, rhMean_results, aridityMean_results, elevation_results)
variables = c('tmin', 'tmax', 'trange', 'precipTot', 'rhMean', 'aridityMean', 'elevation')

makeGEMMAPlots = lapply(1:length(variables), function(i) {
  png(here('Plots', sprintf('GEMMA_univariate_%s_manhattan.png', variables[i])), height = 480, width = 720)
  manhattan(gemma_results[i][[1]], 
            snp = 'X',
            bp = 'BP',
            p = variables[i],
            chr = 'Chr',
            main = sprintf('GEMMA univariate GEA for %s, maf > 0.01', variables[i]),
            highlight = c(inv4mSNPs, hsftf9SNPs))
  dev.off()
})

readGEMMAPlots <- lapply (variables, function(i) {
  rasterGrob(readPNG(file.path(here('Plots', sprintf('GEMMA_univariate_%s_manhattan.png', i))), native = FALSE),
             interpolate = FALSE)
})

pdf(here('Plots', "testgraph.pdf"))
do.call(grid.arrange, c(readGEMMAPlots, ncol = 2, nrow = 4))
dev.off()