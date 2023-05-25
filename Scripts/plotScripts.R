####################################################################################
# Generate plots
#
# Author: Forrest Li
# Plots from analyses
####################################################################################

source(here::here('config.R'))

####################################################################################
### create map for elevation and locations
####################################################################################

terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = -50, right = -80, top = 33), 
                                  zoom = 4, maptype = "terrain-background")

ggmap(terrain_kernels)+
  geom_point(data = cimmyt_grow, 
             aes(x = LongNew, y = LatNew, color = elevation),
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

terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = -40, right = -35, top = 33), 
                                  zoom = 4, maptype = "toner-lite", color = 'bw', messaging = F)
ggmap(terrain_kernels)+
  # geom_point(data = trialinfo, 
  #            aes(x = Trial_longitude, y = Trial_latitude, color = Trial_elevation),
  #            pch=16,alpha=0.5,size=2) + 
  geom_point(data = cimmyt_grow,
             aes(x = LongNew, y = LatNew, color = tmin),
             pch=21,alpha=0.75,size=1, show.legend = F) +
  scale_color_continuous(name = "Temperature", low = 'red', high = 'blue') +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 12)) +
  # ggtitle("Accession locations (n = 3515)") +
  xlab("") + 
  ylab("")

####################################################################################
## histograms of climate variables
####################################################################################

a = ggplot(cimmyt_grow, aes(x = tmax)) +
  geom_histogram(bins = 50, fill = '#78d6eb') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Maximum temperature')

b = ggplot(cimmyt_grow, aes(x = tmin)) +
  geom_histogram(bins = 50, fill = '#c42127') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Minimum temperature')

c = ggplot(cimmyt_grow, aes(x = trange)) +
  geom_histogram(bins = 50, fill = '#d9a74c') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Temperature range in growing season')

d = ggplot(cimmyt_grow, aes(x = precipTot)) +
  geom_histogram(bins = 50, fill = '#1518d1') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Total precipitation over growing season')

e = ggplot(cimmyt_grow, aes(x = elevation)) +
  geom_histogram(bins = 50, fill = '#914fe8') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Elevation')

f = ggplot(cimmyt_grow, aes(x = rhMean)) +
  geom_histogram(bins = 50, fill = '#76c28e') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Mean relative humidity')

g = ggplot(cimmyt_grow, aes(x = aridityMean)) +
  geom_histogram(bins = 50, fill = '#307546') +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=20)) +
  xlab('') +
  ylab('') +
  ggtitle('Mean aridity')

ggarrange(plotlist = list(a, b, c, d, e, f, g),
          nrow = 7,
          ncol = 1)

####################################################################################
## plot climate variable PCA and correlations between variables
####################################################################################

pca_grow = prcomp(finalMat[5:11], scale = T)
autoplot(pca_grow, colour = 'elevation')

pca_grow_transform = prcomp(finalMat_transform[5:11])
autoplot(pca_grow_transform, colour = 'elevation')

ggpairs(finalMat[4:11])
ggpairs(finalMat_transform[4:11])
pcor(finalMat_transform[4:11])

## compare transformation
plot(finalMat_transform$elevation, finalMat$elevation)
plot(finalMat_transform$tmin, finalMat$tmin)
plot(finalMat_transform$tmax, finalMat$tmax)
plot(finalMat_transform$trange, finalMat$trange)
plot(finalMat_transform$precipTot, finalMat$precipTot)
plot(finalMat_transform$aridityMean, finalMat$aridityMean)
plot(finalMat_transform$rhMean, finalMat$rhMean)


####################################################################################
## this is for linguistics project
####################################################################################
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
### manhattan plots for multivariate GEA
####################################################################################
bonferroni = 1/nrow(gea_results)
manhattan(gea_results, 
          snp = 'SNP',
          bp = 'BP',
          main = 'Multivariate GEA across 7 climate variables, maf > 0.01',
          highlight = top_hits_clumped)

manhattan(gea_results %>% filter(maf > 0.05), 
          snp = 'SNP',
          bp = 'BP',
          main = 'Multivariate GEA across 7 climate variables, maf > 0.05')

manhattan(gea_results %>% filter(CHR == 9 & maf > 0.01), 
          snp = 'SNP',
          bp = 'BP',
          # main = 'HSF9, maf > 0.05',
          highlight = hsftf9SNPs,)

manhattan(gea_results, 
          snp = 'SNP',
          bp = 'BP',
          # genomewideline = -log10(bonferroni),
          highlight = top_hits_clumped)

manhattan(gea_results %>% filter(CHR == 4), 
          snp = 'SNP',
          bp = 'BP',
          # main = 'HSF9, maf > 0.05',
          highlight = inv4mSNPs,)

####################################################################################
## plot GEMMA univariate manhattan plots
####################################################################################
gemma_results = list(tmin_results, tmax_results, trange_results, precipTot_results, rhMean_results, aridityMean_results, elevation_results)
variables = c('tmin', 'tmax', 'trange', 'precipTot', 'rhMean', 'aridityMean', 'elevation')

makeGEMMAPlots = lapply(1:length(variables), function(i) {
  png(here(plot_dir, sprintf('GEMMA_univariate_%s_manhattan.png', variables[i])), height = 480, width = 720)
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
  rasterGrob(readPNG(file.path(here(plot_dir, sprintf('GEMMA_univariate_%s_manhattan.png', i))), native = FALSE),
             interpolate = FALSE)
})

pdf(here(plot_dir, "testgraph.pdf"))
do.call(grid.arrange, c(readGEMMAPlots, ncol = 2, nrow = 4))
dev.off()


####################################################################################
## phenotypic consquence of top GEA SNPs
####################################################################################

## general plotting
significance_plot_list = purrr::map2(significance_pvals,
                                     list('bare cob weight', 'grain weight per hectare', 'field weight',
                                          'plant height', 'days to flowering', 'ASI'),
                                     plotpValueSDHorizontal)

significance_fig = ggarrange(plotlist = significance_plot_list, 
                             nrow = 1, ncol = 6,
                             common.legend = F)
annotate_figure(significance_fig, 
                left = textGrob('-log10(p-value)'),
                top = text_grob('SNP significance for phenotype across field trials', face = 'bold', size = 16),
                fig.lab.size = 8)

## specific scripts for talk
significance_pvals = purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process),
                                 list('bare cob weight', 'grain weight per hectare', 'field weight',
                                      'plant height', 'days to flowering', 'ASI'),
                                 plotVariableSignificance)

talk_plot_list = purrr::map2(significance_pvals,
                             list('bare cob weight', 'grain weight per hectare', 'field weight',
                                  'plant height', 'days to flowering', 'ASI'),
                             plotpValueSDTalk)

talk_significance_fig = ggarrange(plotlist = talk_plot_list, 
                                  nrow = 1, ncol = 6,
                                  common.legend = T)

annotate_figure(talk_significance_fig, 
                left = text_grob('Mean -log10(p-value)', size = 14),
                fig.lab.size = 24)

## for poster and figures
poster_plot_list = purrr::map2(list(gwph_process, fw_process, dtf_process, asi_process),
                               list('grain weight per hectare', 'field weight',
                                    'days to flowering', 'ASI'), 
                               plotVariableSignificance)
poster_fig = ggarrange(plotlist = poster_plot_list, 
                       nrow = 4, ncol = 1,
                       common.legend = T)
(fig4 = annotate_figure(poster_fig, 
                        bottom = text_grob('-log10 p-value', size = 16),
                        top = text_grob('SNP significance for phenotype across field trials', size = 18),
                        fig.lab.size = 18,
))

png(here(plot_dir, 'phenotypic-consequence', 'top-SNPs-vs-random.png'), width = 655, height = 655)
print(fig4)
dev.off()