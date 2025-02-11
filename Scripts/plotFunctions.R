####################################################################################
# Plotting functions for SeeDs GEA analysis
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################

library(ggplot2)
library(qqman)
library(ggmap)
library(ggpubr)
library(png)
library(grid)
library(gridExtra)
library(GGally)
library(ggfortify)
library(gridExtra)
library(cowplot)
library(ggtext)
library(lemon)
library(tiff)



############################################################################
## for map plotting
############################################################################
trial_info = read.csv('Phenotype_data/Trial_info.csv')
trial_worldclim = read.csv('Phenotype_data/Trial_worldclim.csv')
trial_master = merge(trial_worldclim, trial_info, by = 'Experimento')

## function that plots manhattan plot
## highlight_SNP_list: pass a list of SNPs to highlight,
## chr_filter: specific chromosomes to plot
## sig: choose threshold for significance line 
## chr_col, bp_col, p_col: column names corresponding to chromosome, BP, and P-value
manual_manhattan = function(results, highlight_SNP_list, chr_filter = NA, sig = 1e-5,
                            chr_col = 'CHR', bp_col = 'BP', p_col = 'P', 
                            title = 'Multivariate GEA'){
  if(!is.na(chr_filter)) {
    results = results %>% filter((!!sym(chr_col)) == chr_filter)
  }
  
  max_BP <- results |>
    group_by((!!sym(chr_col))) |>
    summarise(max_bp = max((!!sym(bp_col)))) |>
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
    dplyr::select((!!sym(chr_col)), bp_add)
  
  results_add <- results |>
    inner_join(max_BP, by = chr_col) |>
    mutate(bp_cum = (!!sym(bp_col)) + bp_add)
  
  axis_set <- results_add |>
    group_by((!!sym(chr_col))) |>
    summarize(center = mean(bp_cum))
  
  ylim <- results_add |>
    filter((!!sym(p_col)) == min((!!sym(p_col)))) |>
    mutate(ylim = abs(floor(log10((!!sym(p_col))))) + 2) |>
    pull(ylim)
  
  # sig <- 1e-5
  
  (manhplot <- ggplot(results_add, aes(
    x = bp_cum, y = -log10((!!sym(p_col))),
    color = as.factor((!!sym(chr_col))), size = -log10((!!sym(p_col)))
  )) +
      geom_hline(
        yintercept = -log10(sig), color = "blue",
        linetype = "dashed"
      ) +
      geom_point(alpha = 0.75, size = 0.5) +
      scale_x_continuous(
        label = axis_set %>% pull((!!sym(chr_col))),
        breaks = axis_set$center
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
      scale_color_manual(values = rep(
        c("grey4", "grey30"),
        unique(length(axis_set %>% pull(!!sym(chr_col))))
      )) +
      geom_point(data = results_add[results_add$SNP %in% highlight_SNP_list,],
                 alpha = 0.75, size = 0.75, color = 'red') + 
      scale_size_continuous(range = c(0.5, 3)) +
      labs(
        x = NULL,
        y = "-log<sub>10</sub>(p)"
      ) +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5)
      ) +
      theme_bw() +
      ggtitle(title)
  )
  
}

### code to easily find SNP names within coordinate intervals for inv4m and hsftf9 for plotting 
{
  tmin_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmin', 'p_wald')
  inv4mCoords = c(1.71e+8, 1.86e8)
  hsftf9Coords = c(1.480e8 , 1.490e8)
  inv4mSNPs = tmin_results %>% filter(BP > inv4mCoords[1] & BP < inv4mCoords[2] & Chr == 4) %>% pull(V1)
  hsftf9SNPs = tmin_results %>% filter(BP > hsftf9Coords[1] & BP < hsftf9Coords[2] & Chr == 9) %>% pull(V1)
}

## returns min, mean, and limits of data for plotting
min.mean.sd.max <- function(x) {
  mean_val = mean(x, na.rm = T)
  r <- c(min(x, na.rm = T), mean_val - sd(x, na.rm = T), mean_val, mean_val + sd(x, na.rm = T), max(x, na.rm = T))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black", size = 0.5) +
    geom_smooth(method = method, color = "red", ...)
  p
}

plotpValueSDFigure = function(df, phenotype) {
  ggplot(df %>% 
           filter(is_top_hit | matching), 
         aes(x = is_top_hit, y = logp, color = is_top_hit, size = is_top_hit)) +
    # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
    geom_jitter(width = 0.15, show.legend = F) +
    ggtitle(str_wrap(phenotype))+
    scale_x_discrete(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs')) +
    scale_color_manual(values = c("#8796c7", "#eb4034")) +
    scale_size_manual(values = c(1, 2)) +
    theme(axis.text = element_text(size = 24),
          axis.ticks = element_text(size = 18),
          legend.position = 'none') +
    theme_bw() +
    coord_flip() +
    ylim(0.3, 0.75) +
    labs(x = '', y = '') +
    theme(plot.title = element_text(size = 16))
}

plotpValueSDHorizontal = function(df, phenotype) {
  ggplot(df %>% 
           filter(is_top_hit | matching), 
         aes(x = is_top_hit, y = logp, color = is_top_hit, size = is_top_hit)) +
    # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
    geom_jitter(width = 0.15, show.legend = F) +
    ggtitle(str_wrap(phenotype))+
    scale_x_discrete(labels = c('bootstrapped\nruns of\n random SNPs', 'top GEA \nSNPs')) +
    scale_color_manual(values = c("#8796c7", "#eb4034")) +
    scale_size_manual(values = c(1, 3)) +
    theme(axis.text = element_text(size = 24),
          axis.ticks = element_text(size = 18),
          legend.position = 'none') +
    theme_bw() +
    coord_flip() +
    ylim(0.25, 0.75) +
    labs(x = '', y = '') +
    theme(plot.title = element_text(size = 14))
}

plotpValueSDTalk = function(df, phenotype) {
  ggplot(df %>% 
           filter(is_top_hit | matching), 
         aes(x = is_top_hit, y = logp, color = is_top_hit, size = is_top_hit)) +
    # stat_summary(fun.data = min.mean.sd.max, geom= 'boxplot')+
    geom_jitter(width = 0.15) +
    ggtitle(str_wrap(phenotype))+
    scale_x_discrete(labels = c('', '')) +
    scale_color_manual(labels = c('bootstrapped runs of random SNPs', 'top GEA SNPs'), values = c("#8796c7", "#eb4034"), ) +
    scale_size_manual(values = c(2, 4)) +
    theme(legend.text = element_text(size = 30),
          axis.text = element_text(size = 30),
          axis.ticks = element_text(size = 30),
          legend.key.size = unit(2, 'cm'),
          plot.title = element_text(size = 30)) +
    guides(scale = 'none', size = 'none') +
    theme_bw() +
    ylim(0.25, 0.75) +
    labs(x = '', y = '', color = '')
}

plotpValueSDManuscript = function(df) {
  ggplot(df %>% 
           filter(is_top_hit | matching) %>% arrange(is_top_hit), 
         aes(x = trait, y = logP, color = is_top_hit, size = is_top_hit)) +
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

## function to obtains SE error bars for plotting
getSEs = function(y) {
  se = sd(y)/sqrt(length(y))
  mu = mean(y)
  c(ymin = mu-se, ymax = mu+se)
}

## shifts legend to be placed in empty panel in multi-facet plot
shift_legend <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

## label names and colors for different models
all_columns = c("lm_PC5", 'lm_matching_SNPs', 'lm_enriched_SNPs', 
                'lm_all_SNPs', 'lm_all_SNPs_env',
                'lm_PC5_env', 'lm_env', 'rf_env', 'rf_env_PC5',
                'lm_unstructured_SNPs_only', 'lm_unstructured_SNPs', 'lm_unstructured_matching_SNPs_only', 'lm_unstructured_matching_SNPs')

all_labels = c('Top 5 PCs', 'Random SNPs + PCs', 'Top envGWAS SNPs + PCs', 
               'All SNPs', 'All SNPs + Climate',
               'PCs + Climate', 'Climate data LMM', 'Climate data RF', 'Climate data + PCs RF', 
               'Uncorrected SNPs', 'Uncorrected SNPs + PCs', 'Matching uncorrected SNPs', 'Matching uncorrected SNPs + PCs')

model_palette = c('#4F4C4D', '#6994dd', '#8db5f1', 
                           '#f9ea71', '#F9C343', 
                           '#E7601F','#F57F47', '#df8b8b', '#CE4646', 
                           '#702963', '#AA336A', '#d02963', '#d00000')

create_prediction_comparison_plot_byModel = function(blup_style, trait_i, columns) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/%s/modelPrediction_nostress_%s_results_resid.csv', blup_style, trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = columns, names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = columns, ordered = TRUE)
  ## shorten gwph name for brevity
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byModel = ggplot(all_results_long, 
                      aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    # geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    # facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
    theme_bw() +
    scale_fill_manual(values = model_palette[which(all_columns %in% columns)], labels = all_labels[which(all_columns %in% columns)]) +
    ylim(c(0, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byModel
}

create_prediction_comparison_plot_byTrial = function(blup_style, trait_i, columns) {
  prediction_accuracy = read.csv(here(sprintf('Analyses/PhenotypicPrediction/%s/modelPrediction_nostress_%s_results_resid.csv', blup_style, trait_i)))
  all_results = merge(prediction_accuracy, trial_master[c('Experimento', 'Localidad', 'Año', 'latitude', 'Trial_elevation', 'meanTemp', 'annualPrecipitation')])
  all_results = all_results[order(all_results$Trial_elevation),]
  all_results$name = paste(all_results$Localidad, all_results$Año)
  all_results$name = factor(all_results$name, levels = all_results$name %>% unique())
  
  ## plot all models together for R predictive accuracy
  all_results_long = pivot_longer(all_results, cols = columns, names_to = 'model')
  all_results_long$model = factor(all_results_long$model, 
                                  levels = columns, ordered = TRUE)
  ## shorten gwph name for brevity
  if(trait_i == 'GrainWeightPerHectareCorrected')
    trait_i = 'GrainWeight'
  plot_byModel = ggplot(all_results_long, 
                        aes(x = model, y = value, fill = model)) +
    stat_summary(fun = mean, geom = 'bar') +
    stat_summary(fun.data = getSEs, geom = 'errorbar', width = 0.1) +
    geom_jitter(aes(x = model, y = value), size = 0.1, height = 0, width = 0.1) +
    facet_wrap(facets = 'Experimento') +
    ylab('Pearson correlation (r)') +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = model_palette[which(all_columns %in% columns)], labels = all_labels[which(all_columns %in% columns)]) +
    ylim(c(0, 0.4)) +
    ggtitle(sprintf('%s', gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", trait_i)))
  plot_byModel
}

rotation_hjust = function(angle) {
  rads = (angle - 0) * pi / 180
  hjust = 0.5 * (1 - sin(rads))
  hjust
}

rotation_vjust = function(angle) {
  rads = (angle - 90) * pi / 180
  vjust = 0.5 * (1 + cos(rads))
  vjust
}
