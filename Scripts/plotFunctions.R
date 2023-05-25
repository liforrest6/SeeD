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

tmin_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmin', 'p_wald')
inv4mCoords = c(1.71e+8, 1.86e8)
hsftf9Coords = c(1.480e8 , 1.490e8)
inv4mSNPs = tmin_results %>% filter(BP > inv4mCoords[1] & BP < inv4mCoords[2] & Chr == 4) %>% pull(X)
hsftf9SNPs = tmin_results %>% filter(BP > hsftf9Coords[1] & BP < hsftf9Coords[2] & Chr == 9) %>% pull(X)

min.mean.sd.max <- function(x) {
  mean_val = mean(x, na.rm = T)
  r <- c(min(x, na.rm = T), mean_val - sd(x, na.rm = T), mean_val, mean_val + sd(x, na.rm = T), max(x, na.rm = T))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

sampleSNPs = function(df, top_hits, number_bootstraps, gea = F) {
  ### function that takes a Joint result from processJointresults and the top hits from envGWAS
  ### and adds columns that identify top SNP hits as well as samples SNPs with matching maf
  # create factor column from maf using cut (20 intervals in format [0, 0.025))
  if(gea == F) {
    df$maf_factor = cut(df$maf, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
  } else if(gea == T) {
    df$maf_factor = cut(df$maf.gea, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
  }
  df$is_top_hit = df$SNP %in% top_hits
  
  runs = data.frame(SNP = df$SNP)
  
  nested_df = df %>% 
    # nest by maf interval factor
    group_by(maf_factor) %>% 
    nest() %>% 
    ungroup() %>%
    # create new column of just top hits
    mutate(top_hits = map(data, ~filter(., is_top_hit))) %>% 
    # calculate frequency of top hits
    mutate(count = map_dbl(top_hits, nrow)) %>% 
    # mutate(freq = round(count / sum(count), 3)) %>% 
    # calculate new column of no top hits
    mutate(non_top_hits = map(data, ~filter(., !is_top_hit)))
  
  for(i in 1:number_bootstraps){
    bootstrap_SNPs = nested_df %>% 
      # create new column by sampling non-top-hits using frequency of top hits
      mutate(samp = map2(non_top_hits, count, sample_n)) %>% 
      dplyr::select(samp) %>% 
      unnest(samp) %>% 
      pull(SNP)
    runs[[sprintf('bootstrap%d', i)]] = F
    runs[[sprintf('bootstrap%d', i)]][which(runs$SNP %in% bootstrap_SNPs)] = T
  }
  return(runs)
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

## combine joint GWAS results to effect sizes, get -log10 on p-values
plotVariableSignificance = function(phenotype, phenotype_name){
  number_bootstraps = 30
  phenotype_two_fx = merge_JGWAS_GEA(phenotype, gea_beta_hats)
  SNP_matrix = sampleSNPs(phenotype_two_fx, top_hits_clumped, number_bootstraps, gea = T)
  # pcollection = data.frame(matrix(ncol = 3, nrow = 0))
  # colnames(pcollection) = c('name', 'is_top_hit', 'p-value')
  pcollection = data.frame('name' = 'gea_hits',
                           'is_top_hit' = T,
                           'matching' = F,
                           'logp' = phenotype_two_fx %>% 
                             filter(SNP %in% top_hits_clumped) %>% 
                             pull(logp) %>% 
                             mean(., na.rm = T))
  for(j in 1:number_bootstraps){
    pvalue = mean(phenotype_two_fx[which(SNP_matrix[[sprintf('bootstrap%d',j)]] == T),]$logp)
    pcollection = rbind(pcollection, c(sprintf('bootstrap%d',j), F, T, pvalue))
  }
  
  pcollection$is_top_hit = as.logical(pcollection$is_top_hit)
  pcollection$matching = as.logical(pcollection$matching)
  pcollection$logp = as.numeric(pcollection$logp)
  # x = plotpValueSDFigure(pcollection, phenotype_name)
  return(pcollection)
}

