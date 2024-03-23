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
inv4mSNPs = tmin_results %>% filter(BP > inv4mCoords[1] & BP < inv4mCoords[2] & Chr == 4) %>% pull(V1)
hsftf9SNPs = tmin_results %>% filter(BP > hsftf9Coords[1] & BP < hsftf9Coords[2] & Chr == 9) %>% pull(V1)

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

## combine joint GWAS results to effect sizes, get -log10 on p-values - do this for fisher p-values
calculateVariableSignificanceFisher = function(phenotype, phenotype_name){
  set.seed(42)
  number_bootstraps = 1000
  phenotype_two_fx = merge_JGWAS_GEA(phenotype, gea_beta_hats)
  SNP_matrix = sampleSNPs(phenotype_two_fx, top_hits_clumped, number_bootstraps, gea = T)
  # pcollection = data.frame(matrix(ncol = 3, nrow = 0))
  # colnames(pcollection) = c('name', 'is_top_hit', 'p-value')
  pcollection = data.frame('name' = 'gea_hits',
                           'is_top_hit' = T,
                           'matching' = F,
                           'fisherP' = phenotype_two_fx %>% 
                             filter(SNP %in% top_hits_clumped) %>% 
                             pull(fisherP) %>% 
                             mean(., na.rm = T))
  
  sampled_pcollection = data.frame(do.call(rbind, lapply(c(1:number_bootstraps), function(j){
    pvalue = mean(phenotype_two_fx[which(SNP_matrix[[sprintf('bootstrap%d',j)]] == T),]$fisherP, na.rm = T)
    this_pcollection = c(sprintf('bootstrap%d',j), F, T, pvalue)
    this_pcollection
  })))
  colnames(sampled_pcollection) = c('name', 'is_top_hit', 'matching', 'fisherP')
  pcollection = rbind(pcollection, data.frame(sampled_pcollection))
  
  
  # for(j in 1:number_bootstraps){
  #   pvalue = mean(phenotype_two_fx[which(SNP_matrix[[sprintf('bootstrap%d',j)]] == T),]$fisherP)
  #   pcollection = rbind(pcollection, c(sprintf('bootstrap%d',j), F, T, pvalue))
  # }
  
  pcollection$is_top_hit = as.logical(pcollection$is_top_hit)
  pcollection$matching = as.logical(pcollection$matching)
  pcollection$fisherP = as.numeric(pcollection$fisherP)
  pcollection$logP = -log10(pcollection$fisherP)
  # x = plotpValueSDFigure(pcollection, phenotype_name)
  return(pcollection)
}

## original calculateVariableSignificance function without Fisher's P value
calculateVariableSignificance = function(phenotype, phenotype_name){
  set.seed(42)
  number_bootstraps = 1000
  phenotype_two_fx = merge_JGWAS_GEA(phenotype, gea_beta_hats)
  SNP_matrix = sampleSNPs(phenotype_two_fx, top_hits_clumped, number_bootstraps, gea = T)
  # pcollection = data.frame(matrix(ncol = 3, nrow = 0))
  # colnames(pcollection) = c('name', 'is_top_hit', 'p-value')
  pcollection = data.frame('name' = 'gea_hits',
                           'is_top_hit' = T,
                           'matching' = F,
                           'logP' = phenotype_two_fx %>% 
                             filter(SNP %in% top_hits_clumped) %>% 
                             pull(logP) %>% 
                             mean(., na.rm = T))
  for(j in 1:number_bootstraps){
    pvalue = mean(phenotype_two_fx[which(SNP_matrix[[sprintf('bootstrap%d',j)]] == T),]$logP)
    pcollection = rbind(pcollection, c(sprintf('bootstrap%d',j), F, T, pvalue))
  }
  
  pcollection$is_top_hit = as.logical(pcollection$is_top_hit)
  pcollection$matching = as.logical(pcollection$matching)
  pcollection$logP = as.numeric(pcollection$logP)
  # x = plotpValueSDFigure(pcollection, phenotype_name)
  return(pcollection)
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
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5)
      ) +
      ggtitle(title)
  )
  
}


sampleSNPs_maf = function(df, top_hits, number_bootstraps, gea = F) {
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
    runs[[sprintf('bootstrap%d', i)]][which(runs$SNP %in% bootstrap_SNPs)] = df[which(runs$SNP %in% bootstrap_SNPs), 'maf_factor']
  }
  return(runs)
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
