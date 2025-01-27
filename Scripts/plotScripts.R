####################################################################################
# Generate plots
#
# Plots from analyses
####################################################################################

source(here::here('config.R'))

####################################################################################
### create map for elevation and locations
####################################################################################


terrain_kernels = get_stadiamap(c(left = -120, bottom = -40, right = -35, top = 33),
                                zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

cimmyt_grow = read.csv(here::here(env_data_dir, 'GEA-climate-nontransformed.csv'))
trial_info = read.csv('Phenotype_data/Trial_info.csv')

(elevation_map = ggmap(terrain_kernels, extent = 'device')+
  geom_point(data = cimmyt_grow, 
             aes(x = LongNew, y = LatNew, color = elevation),
             pch=16,alpha=0.25,size=.5) + 
  scale_color_continuous(name = "Elevation (m)", low = 'blue', high = 'red') +
  geom_point(data = trial_info, aes(x = Trial_longitude, y =Trial_latitude),
             pch = 24, alpha = 1, size = .3, color = '#3da834') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 10)) +
  theme(legend.title = element_text(size = 5),
        legend.text = element_text(size = 4),
        legend.position = c(0.2, 0.2),
        legend.key.size = unit(.2, "lines")) +
  # guides(shape = guide_legend(override.aes = list(size = 6)),
  #        color = guide_legend(override.aes = list(size = 6))) +
  # ggtitle("Geo-coordinates for CIMMYT data", ) +
  xlab("") +
  ylab("")
)

## run transfer_functions.R first before running this
# (environmental_map_transfer_plots = plot_grid(
#   plot_grid(elevation_map, transfer_plots, labels = 'AUTO', nrow = 2, rel_heights = c(2,1))))
# png(here(plot_dir, 'Manuscript', 'environmental_map_transfer_plots.png'), width = 600, height = 600)
# print(environmental_map_transfer_plots)
# dev.off()

pdf(here(plot_dir, 'Manuscript', 'environmental_map.pdf'), width = 3, height = 3)
# ggsave(here(plot_dir, 'Manuscript', 'environmental_map.png'), width = 500, height = 500, dpi = 300)
print(elevation_map)
dev.off()

tiff(here(plot_dir, 'Manuscript', 'environmental_map.tif'), 
     width = 3, height = 3, units = 'in', res = 1000)
print(elevation_map)
# ggsave(here(plot_dir, 'Manuscript', 'environmental_map.png'), width = 500, height = 500, dpi = 300)
dev.off()

ggsave('environmental_map.pdf',
       plot = elevation_map,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 80,
       height = 80,
       units = 'mm'
)

# png(here(plot_dir, 'environmental_map.png'), width = 500, height = 500)
# print(elevation_map)
# dev.off()


(pca_map = ggmap(terrain_kernels)+
    geom_point(data = cimmyt_grow, 
               aes(x = LongNew, y = LatNew, color = elevation),
               pch=16,alpha=0.2,size=1) + 
    scale_color_continuous(name = "Elevation (m)", low = 'darkgreen', high = 'red') +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.position = c(0.15, 0.15)) +
    # guides(shape = guide_legend(override.aes = list(size = 6)),
    #        color = guide_legend(override.aes = list(size = 6))) +
    # ggtitle("Geo-coordinates for CIMMYT data", ) +
    xlab("") +
    ylab("")
)

mexico_terrain_kernels = get_stadiamap(c(left = -120, bottom = 0, right = -75, top = 33),
                                zoom = 4, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

ggmap(mexico_terrain_kernels)+
  geom_point(data = trial_info,
             aes(x = Trial_longitude, y = Trial_latitude, color = Trial_elevation),
             pch=16,alpha=0.5,size=2) +
  # geom_point(data = cimmyt_grow,
  #            aes(x = LongNew, y = LatNew, color = tmin),
  #            pch=21,alpha=0.75,size=1, show.legend = T) +
  scale_color_continuous(name = "min temp", low = 'red', high = 'blue') +
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

(environmental_plots = ggarrange(plotlist = list(a, b, c, d, e, f, g),
          nrow = 7,
          ncol = 1)
)

####################################################################################
## plot climate variable PCA and correlations between variables
####################################################################################
initial_columns = c('Sample.ID', 'LongNew', 'LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
cimmyt_raw = cimmyt_grow %>% dplyr::select(all_of(initial_columns))


(environmental_correlations = ggpairs(cimmyt_raw[4:10], lower = list(continuous = wrap(lowerFn, method = 'loess')),
        upper = list(continuous = wrap('cor', size = 2))) + theme(text = element_text(size = 7))
)

## ggpairs for correlations
ggsave('environmental_correlations.pdf',
       plot = environmental_correlations,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 5,
       height = 5,
       units = 'in'
)

## colored autoplot based on climate variable PCA
## see pca-analysis.R


####################################################################################
### manhattan plots for multivariate GEA
####################################################################################
bonferroni = 1/nrow(gea_results)


manhplot = manual_manhattan(gea_results, top_hits_clumped,title = '')


pdf(here(plot_dir, 'Manuscript', 'gea_qqplot.pdf'), width = 5, height = 5)
(qqman::qq(gea_results %>% filter(maf > 0.05) %>% pull(P), main = 'qqplot, maf > 0.05'))
dev.off()

unstructured_envGWAS = manual_manhattan(unstructured %>% filter(P > 0),
                 highlight_SNP_list = c(), title = '')

manual_manhattan(gea_results,
                 highlight_SNP_list = c(), title = 'structured envGWAS')


qqman::manhattan(unstructured)


####################################################################################
## plot GEMMA univariate manhattan plots
####################################################################################
gemma_results = list(tmin_results, tmax_results, trange_results, precipTot_results, rhMean_results, aridityMean_results, elevation_results)
variables = c('tmin', 'tmax', 'trange', 'precipTot', 'rhMean', 'aridityMean', 'elevation')

# makeGEMMAPlots = lapply(1:length(variables), function(i) {
#   png(here(plot_dir, sprintf('GEMMA_univariate_%s_manhattan.png', variables[i])), height = 480, width = 720)
#   manhattan(gemma_results[i][[1]], 
#             snp = 'X',
#             bp = 'BP',
#             p = variables[i],
#             chr = 'Chr',
#             main = sprintf('GEMMA univariate GEA for %s, maf > 0.01', variables[i]),
#             highlight = c(inv4mSNPs, hsftf9SNPs))
#   dev.off()
# })

## doing this for manual
makeGEMMAPlots = lapply(1:length(variables), function(i) {
  this_gemma = manual_manhattan(gemma_results[i][[1]],
                                bp_col = 'BP',
                                p_col = variables[i],
                                chr_col = 'Chr',
                                title = sprintf('%s GEA', variables[i]),
                                highlight_SNP_list = c(inv4mSNPs, hsftf9SNPs))
  
  png(here(plot_dir, sprintf('GEMMA_univariate_%s_manhattan_manual.png', variables[i])), height = 480, width = 720)
  print(this_gemma)
  dev.off()
  this_gemma
})

readGEMMAPlots <- lapply (variables, function(i) {
  rasterGrob(readPNG(file.path(here(plot_dir, sprintf('GEMMA_univariate_%s_manhattan.png', i))), native = FALSE),
             interpolate = FALSE)
})

# pdf(here(plot_dir, "testgraph.pdf"))
png(here(plot_dir, 'Manuscript', 'univariate_GEA.png'), width = 500, height = 600)
do.call(grid.arrange, c(makeGEMMAPlots, ncol = 2, nrow = 4))
dev.off()


####################################################################################
## enrichment of phenotypic consquence of top GEA SNPs
####################################################################################

## specific scripts for talk
# significance_plot_list = purrr::map2(significance_pvals,
#                                      list('bare cob weight', 'grain weight per hectare', 'field weight',
#                                           'plant height', 'days to flowering', 'ASI'),
#                                      plotpValueSDHorizontal)
# 
# significance_fig = ggarrange(plotlist = significance_plot_list, 
#                              nrow = 1, ncol = 6,
#                              common.legend = F)
# annotate_figure(significance_fig, 
#                 left = textGrob('-log10(p-value)'),
#                 top = text_grob('SNP significance for phenotype across field trials', face = 'bold', size = 16),
#                 fig.lab.size = 8)


# talk_plot_list = purrr::map2(significance_pvals,
#                              list('bare cob weight', 'grain weight per hectare', 'field weight',
#                                   'plant height', 'days to flowering', 'ASI'),
#                              plotpValueSDTalk)
# 
# talk_significance_fig = ggarrange(plotlist = talk_plot_list, 
#                                   nrow = 1, ncol = 6,
#                                   common.legend = T)
# 
# annotate_figure(talk_significance_fig, 
#                 left = text_grob('Mean -log10(p-value)', size = 14),
#                 fig.lab.size = 24)

## plotting GEA enrichment for manuscript
significance_pvals = purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process),
                                 list('bare cob weight', 'grain weight per hectare', 'field weight',
                                      'plant height'),
                                 calculateVariableSignificance)

for(i in 1:4){
  trait_pvals = significance_pvals[[i]]
  lower_pval_count = trait_pvals[trait_pvals$logP >= trait_pvals[1, 'logP'],] %>% count() - 1
  print(lower_pval_count / (nrow(trait_pvals) - 1))
}

# manuscript_plot_list = purrr::map2(significance_pvals,
#                              list('bare cob weight', 'grain weight per hectare', 'field weight',
#                                   'plant height', 'days to flowering', 'ASI'),
#                              plotpValueSDHorizontal)

allTraits_pvals = bind_rows(significance_pvals, .id = 'id')
# allTraits_pvals$trait <- dplyr::recode(allTraits_pvals$id, 
#                                from=c("1","2","3","4","5","6"), 
#                                to=c('bare cob weight',  
#                                     , 'ASI'))
allTraits_pvals = allTraits_pvals %>% mutate(trait = recode(id,
                                           "1" = 'bare cob weight', 
                                           "2" = 'grain weight per hectare',
                                           "3" = 'field weight',
                                           "4" = 'plant height'))

(manuscript_significance_fig = plotpValueSDManuscript(allTraits_pvals)
)

## deprecated for plotting different traits in different facets
# (manuscript_significance_fig = ggarrange(plotlist = manuscript_plot_list, 
#                                   nrow = 6, ncol = 1,
#                                   common.legend = T)
# )


# ## for poster and figures
# poster_plot_list = purrr::map2(list(gwph_process, fw_process, dtf_process, asi_process),
#                                list('grain weight per hectare', 'field weight',
#                                     'days to flowering', 'ASI'), 
#                                plotVariableSignificance)
# poster_fig = ggarrange(plotlist = poster_plot_list, 
#                        nrow = 4, ncol = 1,
#                        common.legend = T)
# (fig4 = annotate_figure(poster_fig, 
#                         bottom = text_grob('-log10 p-value', size = 16),
#                         top = text_grob('SNP significance for phenotype across field trials', size = 18),
#                         fig.lab.size = 18,
# ))

# png(here(plot_dir, 'phenotypic-consequence', 'top-SNPs-vs-random.png'), width = 655, height = 655)
# print(fig4)
# dev.off()
environmental_palette = c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51')

(enrichment_by_violin = ggplot(allTraits_pvals %>% 
         filter(matching), 
       aes(x = trait, y = logP, fill = trait)) +
  # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
  geom_violin(draw_quantiles = T, trim = F) +
  geom_point(data = allTraits_pvals %>% filter(is_top_hit), color = 'black', fill = 'red', size = 2, shape = 23) +
  # ggtitle(str_wrap(phenotype))+
  # scale_x_discrete(angle = 20, size = 12, vjust = 0.5) +
  # scale_color_manual(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs'), values = c("#8796c7", "#eb4034"), ) +
  # scale_size_manual(values = c(2, 4)) +
  theme_bw() +
  labs(x = '', y = '-log(p-value)', color = '') +
  theme(legend.text = element_text(size = 15),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        # axis.ticks = element_text(size = 30),
        # axis.ticks.x = element_text(angle = 20, size = 30, vjust = 0.5),
        # legend.key.size = unit(10, 'cm'),
        plot.title = element_text(size = 20)) +
  scale_fill_manual(values = environmental_palette, guide = 'none') +
  # guides(scale = 'none', size = 'none', color = guide_legend(override.aes = list(size = 6))) +
  # guides(scale = 'none', size = 'none', fill = 'none', shape = c('GEA SNPs')) +
  ylim(0, 0.8) +
  coord_flip()
)

# (gea_and_enrichment = plot_grid(manhplot, plot_grid(enrichment_by_violin, plot_bySNPs, ncol = 2, rel_widths = c(2, 3)), 
#                                 labels = 'AUTO', nrow = 2, rel_heights = c(2,3)))
(gea_and_enrichment = plot_grid(manhplot + theme(axis.title.y = element_markdown(size = 7), axis.text.x = element_text(size = 7)), 
                                enrichment_by_violin, labels = 'AUTO', nrow = 2, rel_widths = c(3, 2)))


ggsave('gea_and_enrichment.pdf',
       plot = gea_and_enrichment,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 80,
       height = 80,
       units = 'mm'
)

####################################################################################
## plot accuracy prediction of models
####################################################################################
## run gea-analysis.R first

trial_worldclim = read.csv('Phenotype_data/Trial_worldclim.csv')
trial_master = merge(trial_worldclim, trial_info, by = 'Experimento')

## plot correlations between trials for all raw worldclim variables
cor_between_trials = data.frame(cor(t(trial_worldclim[4:22])), 
                                trial_worldclim$Experimento)
colnames(cor_between_trials) = c(trial_worldclim$Experimento, 'Trial')

trial_correlations = ggplot(melt(cor_between_trials) %>% arrange(desc(Trial), desc(variable)), 
       aes(x = factor(Trial, level = trial_worldclim$Experimento), y = factor(variable, level = trial_worldclim$Experimento), fill = value)) + 
  geom_tile() +
  scale_fill_gradient2() +
  xlab('') +
  ylab('') +
  labs(fill = 'r2') +
  ggtitle('Correlation between field trial locations') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png(here(plot_dir, 'Manuscript', 'trial_corr.png'), width = 491, height = 491)
print(trial_correlations)
dev.off()



## labels for plotting
traits = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","BareCobWeight","PlantHeight")
trait = traits[3]

## select blups_std or blups_deregressed - we choose deregressed to throw out trials with too much noise
blup_style = 'unstructured_deregressed_nostress_filterTrial'
if(blup_style == 'std'){
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/fiveModels/modelPrediction_%s_results_resid.csv', trait)))
} else if(blup_style == 'deregressed') {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups/modelPrediction_%s_results_resid.csv', trait)))
} else if(blup_style == 'deregressed_nostress') {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress/modelPrediction_nostress_%s_results_resid.csv', trait)))
} else if(blup_style == 'deregressed_filterTrial') {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_filterTrial/modelPrediction_%s_results_resid.csv', trait)))
} else if(blup_style == 'deregressed_nostress_filterTrial') {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait)))
} else if(blup_style == 'unstructured_deregressed_nostress_filterTrial') {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/unstructured_deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait)))
}


## compile all blup prediction results
all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
all_results = all_results[order(all_results$Trial_elevation),]
all_results$name = paste(all_results$Localidad, all_results$Año)
all_results$name = factor(all_results$name, levels = all_results$name %>% unique())

## plot all models together for R predictive accuracy
all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, rf_env_PC5, lm_PC5_env, 
                                                      lm_all_SNPs, lm_all_SNPs_env, lm_env,
                                                      lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
all_results_long$model = factor(all_results_long$model, 
                                levels = c("lm_PC5", 'lm_matching_SNPs', 'lm_env',
                                            'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                            'lm_PC5_env', 'rf_env', 'rf_env_PC5'), ordered = TRUE)



# violin plot
# plot_byAll = ggplot(all_results_long, aes(x = model, y = value, color = model)) +
#   geom_violin() +
#   facet_wrap(facets = 'name') +
#   ylab('Pearson correlation (r)') +
#   theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
#   ggtitle(sprintf('%s predictive ability within tester', trait))
# plot_byAll

# bar plot for all models
plot_byAll = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_matching_SNPs', 'lm_enriched_SNPs', 'lm_all_SNPs',
                                                             'lm_all_SNPs', 'lm_all_SNPs_env')), 
                    aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[1:6], labels = all_labels[1:6]) +
  # scale_fill_discrete(labels = all_labels) +
  ggtitle(sprintf('Predictive ability within tester for %s', trait))
plot_byAll
# (plot_byAll = shift_legend(plot_byAll))

## plot improvement of environment model against just PCs
# environment_accuracy = pivot_longer(all_results, cols = c(lm_PC5, rf_env, lm_env, rf_env_PC5), names_to = 'model')
plot_byEnv = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_env', 'rf_env', 'rf_env_PC5')), 
                    aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[c(1, 7, 8, 9)], labels = all_labels[c(1, 7, 8, 9)]) +
  ggtitle(sprintf('Comparison of environmental-only models for %s', trait))
plot_byEnv

## plots all models for manuscript
(phenotypic_prediction = plot_grid(plot_byAll, plot_byEnv, labels = 'AUTO', nrow = 2, rel_heights = c(3, 2)))
png(here(plot_dir, 'Manuscript', 'phenotypic_prediction.png'), width = 600, height = 600)
grid.draw(phenotypic_prediction)
dev.off()

ggsave('phenotypic_prediction.pdf',
       plot = phenotypic_prediction,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 80,
       units = 'mm'
)

## plot aggregate model predictive ability
aggregate_mean_plot = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_matching_SNPs', 'lm_enriched_SNPs', 'lm_all_SNPs',
                                                             'lm_all_SNPs', 'lm_all_SNPs_env')), 
                    aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  # facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[1:6], labels = all_labels[1:6]) +
  # scale_fill_discrete(labels = all_labels) +
  ggtitle(sprintf('Predictive ability within tester for %s', trait))
aggregate_mean_plot

## plot improvement of different models (all SNPs, enriched SNPs, over just using population structure (first 5 PCs))
all_results$enriched_SNPs_improvement = all_results$lm_enriched_SNPs - all_results$lm_PC5
all_results$all_SNPs_improvement = all_results$lm_all_SNPs - all_results$lm_PC5
all_results$all_SNPs_env_improvement = all_results$lm_all_SNPs_env - all_results$lm_PC5
all_results$matching_SNPs_improvement = all_results$lm_matching_SNPs - all_results$lm_PC5
all_results$lm_PC5_env_improvement = all_results$lm_PC5_env - all_results$lm_PC5

## plots improvement of GEA SNPs against random matching SNPs
SNP_improvement_long = pivot_longer(all_results, 
                                    cols = c(all_SNPs_env_improvement, all_SNPs_improvement,
                                             enriched_SNPs_improvement, matching_SNPs_improvement), 
                                    names_to = 'model')

wilcox.test(SNP_improvement_long %>% 
              filter(model == 'all_SNPs_improvement') %>% pull(value),
            SNP_improvement_long %>% 
              filter(model == 'all_SNPs_env_improvement') %>% pull(value),
            paired = T, alternative = 'two.sided')

## plot envGWAS and matching for Figure 3 when comparing prediction of SNPs
plot_bySNPs = ggplot(SNP_improvement_long %>% filter(model %in% c('enriched_SNPs_improvement', 'matching_SNPs_improvement')), aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  scale_size_continuous() +
  facet_wrap(facets = 'Experimento') +
  ylab("Difference in Pearson's r between SNP model - PC model") +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[c(3, 2)], labels = all_labels[c(3, 2)])
  # ggtitle(sprintf('Predictive ability improvement over PC model for %s', trait))
plot_bySNPs

## plot all improvements over PC model
plot_improvement = ggplot(SNP_improvement_long, aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  scale_size_continuous() +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r) of SNP model - PC model') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[c(4, 5, 2, 3)], labels = all_labels[c(4, 5, 2, 3)]) +
  ggtitle(sprintf('Predictive ability improvement over PC model for %s', trait))
plot_improvement
# plot_bySNPs = shift_legend(plot_bySNPs)

# (phenotypic_prediction = plot_grid(plot_byAll, plot_byEnv, labels = 'AUTO', nrow = 2))
# png(here(plot_dir, 'Manuscript', 'phenotypic_prediction_improvement.png'), width = 600, height = 600)
# grid.draw(plot_improvement)
# dev.off()

## this plot is nice for a talk without introducing GEA
(prediction_noGEA = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'rf_env', 'lm_all_SNPs')), 
                          aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  ggtitle(sprintf('Predictive accuracy for %s model prediction within tester', trait)))

### test prediction accuracy across all traits
all_traits_prediction_plot_list = list()
for(trait_i in traits[-c(1)]){
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, lm_matching_SNPs, lm_enriched_SNPs, lm_all_SNPs,
                                                        lm_all_SNPs, lm_all_SNPs_env), names_to = 'model')
  plot_byAll = ggplot(all_results_long, aes(x = model, y = value, color = model)) +
    geom_violin() +
    facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
    scale_fill_discrete(labels = all_labels) +
    ggtitle(sprintf('Predictive accuracy for %s model prediction within tester', trait_i))
  all_traits_prediction_plot_list = append(all_traits_prediction_plot_list, plot_byAll)
  
}

create_prediction_comparison_plot_byTrial = function(trait_i) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, rf_env_PC5, lm_PC5_env, 
                                                        lm_all_SNPs, lm_all_SNPs_env, 
                                                        lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = c("lm_PC5", 'lm_matching_SNPs', 
                                             'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                             'lm_PC5_env', 'rf_env', 'rf_env_PC5'), ordered = TRUE)
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byAll = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_matching_SNPs', 'lm_enriched_SNPs', 'lm_all_SNPs',
                                                               'lm_all_SNPs', 'lm_all_SNPs_env')), aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), 
          axis.title.x = element_blank(), strip.text.x = element_text(size = 5, margin = margin(3, 0, 3, 0, unit = 'pt')))+
    scale_fill_manual(values = model_palette[1:6], labels = all_labels[1:6]) +
    ylim(c(-0.05, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byAll
}

create_prediction_comparison_plot_aggregate = function(trait_i) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, rf_env_PC5, lm_PC5_env, 
                                                        lm_all_SNPs, lm_all_SNPs_env, 
                                                        lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = c("lm_PC5", 'lm_matching_SNPs', 
                                             'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                             'lm_PC5_env', 'rf_env', 'rf_env_PC5'), ordered = TRUE)
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byAll = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_matching_SNPs', 'lm_enriched_SNPs', 'lm_all_SNPs',
                                                               'lm_all_SNPs', 'lm_all_SNPs_env')), aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    # facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), 
          axis.title.x = element_blank(), strip.text.x = element_text(size = 7, margin = margin(3, 0, 3, 0, unit = 'pt')))+
    scale_fill_manual(values = model_palette[1:6], labels = all_labels[1:6]) +
    ylim(c(0, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byAll
}

create_prediction_comparison_plot_environment_byTrial = function(trait_i) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, rf_env_PC5, lm_PC5_env, 
                                                        lm_all_SNPs, lm_all_SNPs_env, lm_env,
                                                        lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = c("lm_PC5", 'lm_matching_SNPs', 'lm_env',
                                             'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                             'lm_PC5_env', 'rf_env', 'rf_env_PC5'), ordered = TRUE)
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byEnv = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_env', 'rf_env', 'rf_env_PC5')), 
                      aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), 
          axis.title.x = element_blank(), strip.text.x = element_text(size = 5, margin = margin(3, 0, 3, 0, unit = 'pt')))+
    scale_fill_manual(values = model_palette[c(1, 7, 8, 9)], labels = all_labels[c(1, 7, 8, 9)]) +
    ylim(c(-0.05, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byEnv
}

create_prediction_comparison_plot_environment_aggregate = function(trait_i) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, rf_env_PC5, lm_PC5_env, 
                                                        lm_all_SNPs, lm_all_SNPs_env, lm_env,
                                                        lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = c("lm_PC5", 'lm_matching_SNPs', 'lm_env',
                                             'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                             'lm_PC5_env', 'rf_env', 'rf_env_PC5'), ordered = TRUE)
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byEnv = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_env', 'rf_env', 'rf_env_PC5')), 
         aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    # facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = model_palette[c(1, 7, 8, 9)], labels = all_labels[c(1, 7, 8, 9)]) +
    ylim(c(0, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byEnv
}

all_traits_prediction_plot_list_byTrial = lapply(traits[-c(1, 2)], create_prediction_comparison_plot_byTrial)
all_traits_prediction_plot_list_aggregate = lapply(traits[-c(1, 2)], create_prediction_comparison_plot_aggregate)
all_traits_prediction_plot_list_environment_byTrial = lapply(traits[-c(1, 2)], create_prediction_comparison_plot_environment_byTrial)
all_traits_prediction_plot_list_environment_aggregate = lapply(traits[-c(1, 2)], create_prediction_comparison_plot_environment_aggregate)

(all_traits_prediction_plots_aggregate = plot_grid(
  plot_grid(plotlist = lapply(all_traits_prediction_plot_list_aggregate, function(p) {p + theme(legend.position = 'none', text = element_text(size = 6))}), 
                                                          nrow = 1, ncol = 5, labels = c('', 'A', 'B', 'C', 'D'),
                                                          legend = get_legend(all_traits_prediction_plot_list_aggregate[[1]] + 
                                                                                theme(legend.position = 'left',text = element_text(size = 6)))),
  plot_grid(plotlist = lapply(all_traits_prediction_plot_list_environment_aggregate, function(p) {p + theme(legend.position = 'none', text = element_text(size = 6))}), 
            nrow = 1, ncol = 5, labels = c('', 'E', 'F', 'G', 'H'),
            legend = get_legend(all_traits_prediction_plot_list_environment_aggregate[[1]] + theme(legend.position = 'left', text = element_text(size = 6)))),
  nrow = 2
))

ggsave('all_traits_prediction_plots_aggregate.pdf',
       plot = all_traits_prediction_plots_aggregate,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 80,
       units = 'mm'
)

(all_traits_prediction_plots = 
  plot_grid(plotlist = lapply(all_traits_prediction_plot_list_byTrial, function(p) {p + theme(legend.position = 'none', text = element_text(size = 8))}), 
            nrow = 3, ncol = 2, labels = c('', 'A', 'B', 'C', 'D'),
            legend = get_legend(all_traits_prediction_plot_list_byTrial[[1]] + theme(legend.position = 'left', text = element_text(size = 8))))
)
  
ggsave('all_traits_prediction_plots.pdf',
       plot = all_traits_prediction_plots,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 200,
       units = 'mm'
)

(all_traits_prediction_plots_byEnv = 
    plot_grid(plotlist = lapply(all_traits_prediction_plot_list_environment_byTrial, function(p) {p + theme(legend.position = 'none', text = element_text(size = 8))}), 
              nrow = 3, ncol = 2, labels = c('', 'A', 'B', 'C', 'D'),
              legend = get_legend(all_traits_prediction_plot_list_environment_byTrial[[1]] + theme(legend.position = 'left', text = element_text(size = 8))))
)

ggsave('all_traits_prediction_plots_byEnv.pdf',
       plot = all_traits_prediction_plots_byEnv,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 200,
       units = 'mm'
)

# for BAPG talk
plot_byAll = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_enriched_SNPs', 'lm_all_SNPs',
                                                             'lm_all_SNPs_env', 'rf_env')), 
                    aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = model_palette[c(1,3,4,5,8)], labels = all_labels[c(1,3,4,5,8)]) +
  # scale_fill_discrete(labels = all_labels) +
  ggtitle(sprintf('Predictive ability within tester for %s', trait))
plot_byAll



## plots for unstructured 

create_prediction_comparison_plot_unstructured = function(trait_i) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/unstructured_deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = c(lm_PC5, lm_all_SNPs, lm_enriched_SNPs,
                                                        lm_enriched_SNPs_only, lm_matching_SNPs,
                                                        lm_unstructured_SNPs, lm_unstructured_SNPs_only), names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = c("lm_PC5", 'lm_matching_SNPs',
                                             'lm_enriched_SNPs', 'lm_all_SNPs',
                                             'lm_unstructured_SNPs_only', 'lm_unstructured_SNPs'), ordered = TRUE)
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byEnv = ggplot(all_results_long %>% filter(model %in% c('lm_PC5', 'lm_enriched_SNPs', 'lm_all_SNPs', 
                                                               'lm_unstructured_SNPs_only', 'lm_unstructured_SNPs')), 
                      aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    # facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = model_palette[c(1, 3, 4, 10, 11)], labels = all_labels[c(1, 3, 4, 10, 11)]) +
    ylim(c(0, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byEnv
}

all_traits_prediction_plot_list_unstructured = lapply(traits[-c(1, 2)], create_prediction_comparison_plot_unstructured)

(all_traits_prediction_plots_unstructured = 
    plot_grid(plotlist = lapply(all_traits_prediction_plot_list_unstructured, function(p) {p + theme(legend.position = 'none', text = element_text(size = 8))}), 
              nrow = 2, ncol = 2, labels = c('B', 'C', 'D', 'E'))
)

unstructured_results = plot_grid(unstructured_envGWAS, 
                                 all_traits_prediction_plots_unstructured, 
                                 get_legend(all_traits_prediction_plot_list_unstructured[[1]] + theme(legend.position = 'bottom', text = element_text(size = 8))), 
                                 nrow = 3, rel_heights = c(2, 3, .5), labels = c('A', ''))

pdf(here(plot_dir, 'Manuscript', 'unstructured_envGWAS.pdf'), width = 6, height = 7)
print(unstructured_results)
dev.off()
# unstructured_long = pivot_longer(all_results, cols = c(lm_PC5, lm_all_SNPs, lm_enriched_SNPs,
#                                                        lm_enriched_SNPs_only, lm_matching_SNPs,
#                                                        lm_unstructured_SNPs, lm_unstructured_SNPs_only), names_to = 'model')

## per trial for unstructured
plot_unstructured = ggplot(unstructured_long, 
                    aes(x = model, y = value, fill = model)) +
  stat_summary(fun = mean, geom = 'bar') +
  stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
  # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
  facet_wrap(facets = 'Experimento') +
  ylab('Pearson correlation (r)') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  # scale_fill_manual(values = model_palette[c(1,3,4,5,8)], labels = all_labels[c(1,3,4,5,8)]) +
  # scale_fill_discrete(labels = all_labels) +
  ggtitle(sprintf('Predictive ability within tester for %s', trait))
plot_unstructured

# read_prediction_accuracy = function(trait_i) {
#   prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/fiveModels/modelPrediction_%s_results_resid.csv', trait_i)))
#   all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
#   all_results = all_results[order(all_results$Trial_elevation),]
#   all_results$name = paste(all_results$Localidad, all_results$Año)
#   all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
#   all_results$trait = trait_i
#   all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, lm_all_SNPs, lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
#   all_results_long = all_results_long %>% mutate(model = recode(model,
#                                                               "lm_PC5" = 'Top 5 PCs', 
#                                                               "rf_env" = 'Enviromic',
#                                                               "lm_all_SNPs" = 'All SNPs',
#                                                               "lm_enriched_SNPs" = 'GEA enriched SNPs only',
#                                                               "lm_matching_SNPs" ='Random SNPs'))
#   
#   all_results_long
# }
# 
# all_accuracies_across_traits = do.call(rbind, lapply(traits[-c(1,5)], read_prediction_accuracy))
# 
# plot_accuracy_by_trait = ggplot(all_accuracies_across_traits, aes(x = trait, y = value, color = trait)) +
#   geom_boxplot() +
#   facet_grid(. ~ model) +
#   ylab('Pearson correlation (r)') +
#   theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())+
#   # scale_color_discrete(labels = all_labels) +
#   ggtitle('Prediction accuracy by trait and model type')
# 
# plot_accuracy_by_trait
# all_traits_prediction_plots = plot_grid(plotlist = lapply(all_traits_prediction_plot_list, function(p) {p + theme(legend.position = 'none')}), 
#                                         nrow = 2, ncol = 3, labels = 'AUTO',
#                                         legend = get_legend(all_traits_prediction_plot_list[[1]] + theme(legend.position = 'left')))
# png(here(plot_dir, 'Manuscript', 'all_traits_prediction_plots.png'), width = 491, height = 982)
# print(all_traits_prediction_plots)
# dev.off()

####################################################################################
## Track top GEA allele geographical distribution
####################################################################################
library(ggmap)
library(leaflet)
library(sp)
finalMat = read.csv(here(env_data_dir, 'GEA-climate-nontransformed.csv'))
GEA_enriched = read.csv(here('Genetic_data/Imputed_V4/GEA_clumped_SNPs_genotype.csv'))
climateAllele = merge(finalMat, GEA_enriched, by.x = 'Unique.ID', by.y = 'X')
climateAllele

# map <- leaflet() %>%
#   addTiles(
#     urlTemplate = "https://tiles.stadiamaps.com/tiles/{variant}/{z}/{x}/{y}{r}.png?api_key=407d34ae-637b-48d8-9c55-87ac9068f78b",
#     attribution = paste('&copy; <a href="https://www.stadiamaps.com/" target="_blank">Stadia Maps</a> ' ,
#                         '&copy; <a href="https://www.stamen.com/" target="_blank">Stamen Design</a> ' ,
#                         '&copy; <a href="https://openmaptiles.org/" target="_blank">OpenMapTiles</a> ' ,
#                         '&copy; <a href="https://www.openstreetmap.org/about" target="_blank">OpenStreetMap</a> contributors'),
#     options = tileOptions(variant='stamen_toner_lite', apikey = 'YOUR-API-KEY')
#   )  %>%
#   fitBounds(lng1 = -86.1581, lat1 = 39.7684, lng2 = -87.1581, lat2 = 40.7684)

df = data.frame(longitude = climateAllele$LongNew,
                latitude = climateAllele$LatNew,
                inv4m = as.factor(climateAllele$S4_177835031),
                hsftf9 = as.factor(climateAllele$S9_148365695),
                chr2locus = as.factor(climateAllele$S2_198537913))
coordinates(df) <- ~longitude+latitude
pal = colorFactor(c('#73BFB8', '#3DA5D9', '#2364AA'), climateAllele$S2_198537913)
leaflet(df) %>% addProviderTiles() %>% addCircleMarkers(color = ~pal(chr2locus), fillOpacity = 0.2, radius = 0.4)

inv4m_map = ggmap(terrain_kernels)+
  # geom_point(data = trialinfo, 
  #            aes(x = Trial_longitude, y = Trial_latitude, color = Trial_elevation),
  #            pch=16,alpha=0.5,size=2) + 
  geom_point(data = climateAllele,
             aes(x = LongNew, y = LatNew, color = as.factor(S4_177835031)),
             pch=21,alpha=0.75,size=1, show.legend = T) +
  scale_color_manual(name = "Inv4m", values = c('blue', 'purple', 'red'), labels = c('Ref', 'Het', 'Alt')) +
  theme(plot.title = element_text(size = 12)) +
  ggtitle("Presence of Inv4m allele") +
  xlab("") + 
  ylab("")

hsftf9_map = ggmap(terrain_kernels)+
  # geom_point(data = trialinfo, 
  #            aes(x = Trial_longitude, y = Trial_latitude, color = Trial_elevation),
  #            pch=16,alpha=0.5,size=2) + 
  geom_point(data = climateAllele,
             aes(x = LongNew, y = LatNew, color = as.factor(S9_148365695)),
             pch=21,alpha=0.75,size=1, show.legend = T) +
  scale_color_manual(name = "Hsftf9", values = c('blue', 'purple', 'red'), labels = c('Ref', 'Het', 'Alt')) +
  theme(plot.title = element_text(size = 12)) +
  ggtitle("Presence of Hsftf9 allele") +
  xlab("") + 
  ylab("")

(allele_spatial_distributions = plot_grid(inv4m_map, hsftf9_map, labels = 'AUTO', nrow = 2))

ggsave('allele_spatial_distributions.pdf',
       plot = allele_spatial_distributions,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 200,
       units = 'mm'
)



####################################################################################
## GEA heatmap of loci vs environmental variable effect sizes
####################################################################################
var_lookup = c(tmin = 'tmin..X', tmax = 'tmax..X', trange = 'trange..X', precipTot = 'precipTot..X', 
           aridityMean = 'aridityMean..X', rhMean = 'rhMean..X', elevation = 'elevation..X')
gea_fx_long = gea_beta_hats %>% 
  rename(all_of(var_lookup)) %>% 
  pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
               names_to = 'env_var')

## plot heatmap of GEA SNP effect sizes, simple version
gea_heatmap = ggplot(gea_fx_long %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = SNP, fill = value)) +
  geom_tile() +
  scale_fill_gradient2() +
  ggtitle('Top GEA-associated SNPs effect size')

gea_gwph_merged = merge_JGWAS_GEA(gwph_process, gea_beta_hats)
gea_gwph_merged$maf_factor = cut(df$maf.gea, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
gea_gwph_merged = gea_gwph_merged %>%  
  rename(all_of(var_lookup)) %>% 
  pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
               names_to = 'env_var')
(GEA_heat_map = ggplot(gea_gwph_merged %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = maf_factor, fill = abs(value))) +
    geom_tile() +
    scale_fill_gradient2(limits = c(0, 0.25)) +
    ggtitle('GEA SNPs beta effect size by maf factor'))

maf_factor_map = data.frame(numbers = c(1:20), factors = gea_gwph_merged$maf_factor %>% unique() %>% sort())


## looking at matching maf SNPs for heat map comparison, grouped by maf factor
matching_SNPs = sampleSNPs_maf(merge_JGWAS_GEA(gwph_process, gea_beta_hats), top_hits_clumped, 50, gea = T)
matching_matrix = matching_SNPs %>% filter(if_any(starts_with('bootstrap'), ~ . != 0))

bootstrap_by_maf = merge(gea_fx_long, matching_matrix, by.x = 'SNP', by.y = 'SNP') %>% 
  pivot_longer(cols = starts_with('bootstrap'),names_to = 'bootstrap_num', values_to = 'maf_factor') %>% 
  mutate(maf_factor = maf_factor_map$factors[match(maf_factor, maf_factor_map$numbers)])

(GEA_heat_map_matching = ggplot(bootstrap_by_maf %>% filter(!is.na(maf_factor)), aes(x = env_var, y = maf_factor, fill = abs(value))) +
    geom_tile() +
    scale_fill_gradient2(limits = c(0, 0.25)) +
    ggtitle('Random matched SNPs beta effect size by maf factor'))

(gea_and_matching_heatmap = plot_grid(GEA_heat_map, GEA_heat_map_matching, labels = 'AUTO', nrow = 2))
png(here(plot_dir, 'Manuscript', 'gea_heatmap.png'), width = 500, height = 750)
print(gea_and_matching_heatmap)
dev.off()
# 
# significance_pvals = purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process),
#                                  list('bare cob weight', 'grain weight per hectare', 'field weight',
#                                       'plant height', 'days to flowering', 'ASI'),
#                                  calculateVariableSignificance)

get_heat_map = function(trait_process, trait_name){
  trait_merged = merge_JGWAS_GEA(trait_process, gea_beta_hats)
  trait_merged$maf_factor = cut(trait_merged$maf.jgwas, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
  trait_merged = trait_merged %>%
    filter(SNP %in% top_hits_clumped) %>% 
    rename(all_of(var_lookup)) %>% 
    pivot_longer(cols = starts_with('Experimentom'), 
                 names_to = 'trial_name',
                 values_to = 'trial_value')
  trait_merged$trait_name = trait_name
  trait_merged
}

all_traits_fx = do.call(rbind, purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process), 
                                           list('bare cob weight', 'grain weight per hectare', 'field weight', 'plant height', 'days to flowering', 'ASI'),
                                           get_heat_map))

(trait_heat_map = ggplot(all_traits_fx %>% filter(SNP %in% top_hits_clumped), aes(x = trait_name, y = maf_factor, fill = abs(trial_value))) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle('GEA SNPs GWPH beta effect size by maf factor'))

ggplot(fw_process %>% pivot_longer(cols = starts_with('Experimentom'),
                                     names_to = 'trial_name',
                                     values_to = 'trial_value'),
       aes(x = trial_name, y = trial_value)) +
  geom_boxplot()

## tried scaling by SD of value for interpretation sake but this doesn't make sense
# gea_sd_scale = gea_beta_hats %>% mutate(tmin..X = tmin..X * sd(finalMat$tmin),
#                                         tmax..X = tmax..X * sd(finalMat$tmax),
#                                         trange..X = trange..X * sd(finalMat$trange),
#                                         aridityMean..X = aridityMean..X * sd(finalMat$aridityMean),
#                                         precipTot..X = precipTot..X * sd(finalMat$precipTot),
#                                         rhMean..X = rhMean..X * sd(finalMat$rhMean),
#                                         elevation..X = elevation..X * sd(finalMat$elevation)
#                                         )
#   
# gea_sd_fx_long = gea_sd_scale %>% 
#   rename(all_of(var_lookup)) %>% 
#   pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
#                names_to = 'env_var')
# (GEA_sd_heat_map = ggplot(gea_sd_fx_long %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = SNP, fill = value)) +
#     geom_tile() +
#     scale_fill_gradient2(limits = c(-0.25, 0.25)) +
#     ggtitle('Top GEA-associated SNPs beta effect size'))
# 
# matching_SNPs = sampleSNPs(merge_JGWAS_GEA(gwph_process, gea_beta_hats), top_hits_clumped, 20, gea = T)
# matching_matrix = matching_SNPs %>% filter(bootstrap1 != 0)
# 
# (GEA_sd_heat_map_matching = ggplot(gea_sd_fx_long %>% filter(SNP %in% matching_matrix$SNP), aes(x = env_var, y = SNP, fill = value)) +
#     geom_tile() +
#     scale_fill_gradient2(limits = c(-0.25, 0.25)) +
#     ggtitle('Random matched SNPs beta effect size'))




## for GSR
ggplot(all_results_long %>% filter(model != 'lm_matching_SNPs', name == 'Agua Fria 2012'), aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_wrap(facets = 'name') +
  ylab('Pearson correlation (r)') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  scale_color_manual(labels = c('Genomic prediction', 'GEA SNP prediction', 'Population structure prediction', 'Enviromic prediction'),
                     values = c('red', 'orange', 'green', 'blue')) +
  ggtitle('Model prediction accuracy for grain weight per hectare')

(transfer_plots = cowplot::plot_grid(plotlist = plots[c(3,5,6)],nrow = 1, ncol=3))
ggarrange(plots[c(3,5,6)], align = 'h', common.legend = T)

ggarrange(plots[c(3,5,6)])
ggarrange(plots[[3]], plots[[5]], plots[[6]], ncol = 3, common.legend = T, legend = 'right')


plotpValueSDManuscriptHorizontal = function(df) {
  ggplot(df %>% 
           filter(matching), 
         aes(x = trait, y = logP)) +
    # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
    geom_violin(draw_quantiles = t) +
    geom_jitter(df %>% filter(is_top_hit), aes(x = trait, y = logP), fill = '#eb4034')
    # ggtitle(str_wrap(phenotype))+
    # scale_x_discrete(angle = 20, size = 12, vjust = 0.5) +
    # scale_color_manual(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs'), values = c("#8796c7", "#eb4034"), ) +
    # scale_size_manual(values = c(2, 4)) +
    theme_bw() +
    labs(x = '', y = '-log(p-value)', color = '') +
    theme(legend.text = element_text(size = 15),
          legend.position = 'top',
          axis.text.x = element_text(angle = 45, size = 15, vjust = 0.5),
          axis.text.y = element_text(size = 15),
          # axis.ticks = element_text(size = 30),
          # axis.ticks.x = element_text(angle = 20, size = 30, vjust = 0.5),
          # legend.key.size = unit(10, 'cm'),
          plot.title = element_text(size = 20)) +
    guides(scale = 'none', size = 'none', color = guide_legend(override.aes = list(size = 6))) +
    ylim(0.25, 0.75) +
    coord_flip()
  
}

plotpValueSDManuscriptHorizontal(allTraits_pvals)

ggplot(allTraits_pvals %>% 
         filter(matching), 
       aes(x = trait, y = logP, fill = trait)) +
  # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
  geom_violin(draw_quantiles = T, trim = F) +
  geom_point(data = allTraits_pvals %>% filter(is_top_hit), color = 'black', size = 5, shape = 9) +
# ggtitle(str_wrap(phenotype))+
# scale_x_discrete(angle = 20, size = 12, vjust = 0.5) +
# scale_color_manual(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs'), values = c("#8796c7", "#eb4034"), ) +
# scale_size_manual(values = c(2, 4)) +
  theme_bw() +
  labs(x = '', y = '-log(p-value)', color = '') +
  theme(legend.text = element_text(size = 15),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, size = 15, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        # axis.ticks = element_text(size = 30),
        # axis.ticks.x = element_text(angle = 20, size = 30, vjust = 0.5),
        # legend.key.size = unit(10, 'cm'),
        plot.title = element_text(size = 20)) +
  # scale_color_discrete(guide = 'none') + 
  # guides(scale = 'none', size = 'none', color = guide_legend(override.aes = list(size = 6))) +
  guides(scale = 'none', size = 'none', fill = 'none', shape = c('GEA SNPs')) +
  ylim(0, 0.8) +
  coord_flip()


