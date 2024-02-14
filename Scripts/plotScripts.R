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
library(ggmap)
library(ggpubr)
library(GGally)
library(gridExtra)
library(cowplot)
library(ggtext)
# library(plyr)
library(dplyr)

stadia_api_key = '407d34ae-637b-48d8-9c55-87ac9068f78b'
register_stadiamaps(stadia_api_key)

terrain_kernels = get_stadiamap(c(left = -120, bottom = -40, right = -35, top = 33),
                                zoom = 5, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

cimmyt_grow = read.csv(here::here(env_data_dir, 'GEA-climate-nontransformed.csv'))

(elevation_map = ggmap(terrain_kernels)+
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
terrain_kernels <- get_stamenmap( bbox = c(left = -120, bottom = -40, right = -35, top = 33), 
                                  zoom = 4, maptype = "toner-lite", color = 'bw', messaging = F)
ggmap(terrain_kernels)+
  geom_point(data = cimmyt_grow,
             aes(x = LongNew, y = LatNew, color = tmin),
             pch=21,alpha=0.75,size=1, show.legend = T) +
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
# ggarrange(plotlist = list(elevation_map, environmental_plots),
#           nrow = 1,
#           ncol = 2)



####################################################################################
## plot climate variable PCA and correlations between variables
####################################################################################
initial_columns = c('Sample.ID', 'LongNew', 'LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
cimmyt_raw = cimmyt_grow %>% dplyr::select(all_of(initial_columns))

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black", size = 0.5) +
    geom_smooth(method = method, color = "red", ...)
  p
}


(environmental_correlations = ggpairs(cimmyt_raw[4:10], lower = list(continuous = wrap(lowerFn, method = 'loess')),
        upper = list(continuous = wrap('cor', size = 4)))
)


pca_grow = prcomp(cimmyt_raw[4:10], scale = T)
autoplot(pca_grow, colour = 'elevation')

plot(pca_grow, type="l")
biplot(pca_grow, xlabs = rep('.', nrow(cimmyt_raw)))

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

## also add PCA plot
(figure1 = plot_grid(
  plot_grid(elevation_map, ggmatrix_gtable(environmental_correlations), labels = 'AUTO', ncol = 2),
  fig5, nrow = 2, rel_heights = c(2,1)))
# figure1 = plot_grid(elevation_map, environmental_correlations, labels = 'AUTO', col = 1, row = 2)
png(here(plot_dir, 'Manuscript', 'environmental_map_transfer_plots.png'), width = 982, height = 982)
print(figure1)
dev.off()

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
          highlight = hsftf9SNPs,
          xlim = c(1.475e8, 1.495e8),
          suggestiveline = F)

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

max_BP <- gea_results %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(CHR, bp_add)

gea_results_add <- gea_results |>
  inner_join(max_BP, by = "CHR") |>
  mutate(bp_cum = BP + bp_add)

axis_set <- gea_results_add |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

ylim <- gea_results_add |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P))) + 2) |>
  pull(ylim)

sig <- 1e-5

(manhplot <- ggplot(gea_results_add, aes(
  x = bp_cum, y = -log10(P),
  color = as.factor(CHR), size = -log10(P)
)) +
  geom_hline(
    yintercept = -log10(sig), color = "blue",
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.75, size = 0.5) +
  scale_x_continuous(
    label = axis_set$CHR,
    breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(
    c("grey4", "grey30"),
    unique(length(axis_set$CHR))
  )) +
  geom_point(data = gea_results_add[gea_results_add$SNP %in% top_hits_clumped,],
             alpha = 0.75, size = 0.75, color = 'red') + 
  scale_size_continuous(range = c(0.5, 3)) +
  labs(
    x = NULL,
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5)
  ) +
  ggtitle('Multivariate GEA')
)

manual_manhattan(gea_results, top_hits_clumped)


(supp1 = qqman::qq(gea_results %>% filter(maf > 0.05) %>% pull(P), main = 'qqplot, maf > 0.05'))
png(here(plot_dir, 'Manuscript', 'supp1-qqplot.png'), width = 491, height = 491)
print(supp1)
dev.off()
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
significance_pvals = purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process),
                                 list('bare cob weight', 'grain weight per hectare', 'field weight',
                                      'plant height', 'days to flowering', 'ASI'),
                                 calculateVariableSignificance)

for(i in 1:6){
  trait_pvals = significance_pvals[[i]]
  lower_pval_count = trait_pvals[trait_pvals$logp >= trait_pvals[1, 'logp'],] %>% count() - 1
  print(lower_pval_count / (nrow(trait_pvals) - 1))
}

manuscript_plot_list = purrr::map2(significance_pvals,
                             list('bare cob weight', 'grain weight per hectare', 'field weight',
                                  'plant height', 'days to flowering', 'ASI'),
                             plotpValueSDHorizontal)

allTraits_pvals = bind_rows(significance_pvals, .id = 'id')
# allTraits_pvals$trait <- dplyr::recode(allTraits_pvals$id, 
#                                from=c("1","2","3","4","5","6"), 
#                                to=c('bare cob weight',  
#                                     , 'ASI'))
allTraits_pvals = allTraits_pvals %>% mutate(trait = recode(id,
                                           "1" = 'bare cob weight', 
                                           "2" = 'grain weight per hectare',
                                           "3" = 'field weight',
                                           "4" = 'plant height',
                                           "5" ='days to flowering',
                                           "6" = 'ASI'))



plotpValueSDManuscript = function(df) {
  ggplot(df %>% 
           filter(is_top_hit | matching) %>% arrange(is_top_hit), 
         aes(x = trait, y = logp, color = is_top_hit, size = is_top_hit)) +
    # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
    geom_jitter(width = 0.15) +
    # ggtitle(str_wrap(phenotype))+
    # scale_x_discrete(angle = 20, size = 12, vjust = 0.5) +
    scale_color_manual(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs'), values = c("#8796c7", "#eb4034"), ) +
    scale_size_manual(values = c(2, 4)) +
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
    ylim(0.25, 0.75)

}

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

(figure2 = plot_grid(plot_grid(manhplot, manuscript_significance_fig, labels = 'AUTO', nrow = 2),
                     GEA_heat_map, ncol = 2, labels = 'AUTO'))
png(here(plot_dir, 'Manuscript', 'gea-and-enrichment.png'), width = 982, height = 982)
print(figure2)
dev.off()

####################################################################################
## plot accuracy prediction of models
####################################################################################
blups = read.csv(here(phenotype_data_dir, 'blups_std.csv'))

accessions_per_trial = read.csv(here(phenotype_data_dir, 'accessions_per_trial_per_trait.csv'))
trial_worldclim = read.csv('Phenotype_data/Trial_worldclim.csv')
trial_info = read.csv('Phenotype_data/Trial_info.csv')
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

png(here(plot_dir, 'Manuscript', 'supp_trial_corr.png'), width = 491, height = 491)
print(trial_correlations)
dev.off()


## summarize blups by tester
blups %>% dplyr::group_by(Trait, Experimento, Tester) %>% dplyr::summarize() %>% count()
count_byTester = blups %>% dplyr::group_by(Trait, Experimento, Tester) %>% dplyr::summarize(n = n())

trait = 'GrainWeightPerHectareCorrected'
prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/fiveModels/modelPrediction_%s_results_resid.csv', trait)))
all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
all_results = all_results[order(all_results$Trial_elevation),]
all_results$name = paste(all_results$Localidad, all_results$Año)
all_results$name = factor(all_results$name, levels = all_results$name %>% unique())

## plot all models together for R predictive accuracy
all_results_long = pivot_longer(all_results, cols = c(lm_PC5, rf_env, lm_all_SNPs, lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
plot_byAll = ggplot(all_results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'name') +
  ylab('Pearson correlation (r)') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  ggtitle(sprintf('Predictive accuracy for %s model prediction within tester', trait))
plot_byAll

## plot improvement of different models (all SNPs, enriched SNPs, over just using population structure (first 5 PCs))
all_results$enriched_SNPs_improvement = all_results$lm_enriched_SNPs - all_results$lm_PC5
all_results$all_SNPs_improvement = all_results$lm_all_SNPs - all_results$lm_PC5
all_results$matching_SNPs_improvement = all_results$lm_matching_SNPs - all_results$lm_PC5

SNP_improvement_long = pivot_longer(all_results, cols = c(all_SNPs_improvement, enriched_SNPs_improvement, matching_SNPs_improvement), names_to = 'model')
plot_bySNPs = ggplot(SNP_improvement_long, aes(x = model, y = value, color = model)) +
  geom_jitter() +
  scale_size_continuous() +
  facet_wrap(facets = 'name') +
  ylab('Pearson correlation (r) of SNP model - PC model') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  ggtitle(sprintf('Improvement on PC model for %s', trait))
plot_bySNPs

## plot improvement of environment model against just PCs
environment_accuracy = pivot_longer(all_results, cols = c(lm_PC5, rf_env), names_to = 'model')
plot_byEnv = ggplot(environment_accuracy, aes(x = model, y = value, color = model, size = n_train)) +
  geom_jitter() +
  facet_wrap(facets = 'name') +
  ylab('Pearson correlation (r)') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  ggtitle(sprintf('Comparing environment model to PC model for %s model', trait))
plot_byEnv

(figure3 = ggarrange(plot_byAll, plot_bySNPs, labels = 'AUTO', ncol = 2))
png(here(plot_dir, 'Manuscript', 'prediction_accuracy.png'), width = 982, height = 491)
print(figure3)
dev.off()

## for tri-lab talk
all_labels = c('All SNPs', 'GEA Enriched SNPs only', 'Random SNPs', 'Top 5 PCs', 'Climate data')
models_for_slide = all_results_long %>% filter(model %in% c('lm_enriched_SNPs', 'lm_matching_SNPs', 'rf_env'))
ggplot(all_results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'name') +
  ylab('Pearson correlation (r)') +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  scale_color_discrete(labels = all_labels) +
  ggtitle('Predictive accuracy for grain weight yield')

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

ggmap(terrain_kernels)+
  # geom_point(data = trialinfo, 
  #            aes(x = Trial_longitude, y = Trial_latitude, color = Trial_elevation),
  #            pch=16,alpha=0.5,size=2) + 
  geom_point(data = climateAllele,
             aes(x = LongNew, y = LatNew, color = S4_177835031),
             pch=21,alpha=0.75,size=1, show.legend = F) +
  scale_color_continuous(name = "hsftf9", low = 'red', high = 'blue') +
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
## GEA heatmap of loci vs environmental variable effect sizes
####################################################################################
var_lookup = c(tmin = 'tmin..X', tmax = 'tmax..X', trange = 'trange..X', precipTot = 'precipTot..X', 
           aridityMean = 'aridityMean..X', rhMean = 'rhMean..X', elevation = 'elevation..X')
gea_fx_long = gea_beta_hats %>% 
  rename(all_of(var_lookup)) %>% 
  pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
               names_to = 'env_var')
gea_fx_long$gv = 2* gea_fx_long$value^2 * gea_fx_long$maf * (1 - gea_fx_long$maf)
(GEA_heat_map = ggplot(gea_fx_long %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = SNP, fill = value)) +
  geom_tile() +
  scale_fill_gradient2() +
  ggtitle('Top GEA-associated SNPs beta effect size'))

GEA_gv_heat_map = ggplot(gea_fx_long %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = SNP, fill = gv)) +
  geom_tile() +
  scale_fill_gradient2() +
  ggtitle('Top GEA-associated SNPs as contribution to genetic variance')

gea_sd_scale = gea_beta_hats %>% mutate(tmin..X = tmin..X * sd(finalMat$tmin),
                                        tmax..X = tmax..X * sd(finalMat$tmax),
                                        trange..X = trange..X * sd(finalMat$trange),
                                        aridityMean..X = aridityMean..X * sd(finalMat$aridityMean),
                                        precipTot..X = precipTot..X * sd(finalMat$precipTot),
                                        rhMean..X = rhMean..X * sd(finalMat$rhMean),
                                        elevation..X = elevation..X * sd(finalMat$elevation)
                                        )
  
gea_sd_fx_long = gea_sd_scale %>% 
  rename(all_of(var_lookup)) %>% 
  pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
               names_to = 'env_var')
(GEA_sd_heat_map = ggplot(gea_sd_fx_long %>% filter(SNP %in% top_hits_clumped), aes(x = env_var, y = SNP, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-0.25, 0.25)) +
    ggtitle('Top GEA-associated SNPs beta effect size'))

matching_SNPs = sampleSNPs(merge_JGWAS_GEA(gwph_process, gea_beta_hats), top_hits_clumped, 1, gea = T)
matching_matrix = matching_SNPs %>% filter(bootstrap1 == T)

(GEA_sd_heat_map_matching = ggplot(gea_sd_fx_long %>% filter(SNP %in% matching_matrix$SNP), aes(x = env_var, y = SNP, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-0.25, 0.25)) +
    ggtitle('Random matched SNPs beta effect size'))
