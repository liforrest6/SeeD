####################################################################################
# Analysis functions for SeeDs GEA analysis
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################

# source(here::here('config.R'))

### compiles files of SNP allele frequencies across all chromosomes
### RETURNS file of all SNPs > 0.05 maf for keeping in GWAS results
# read_MAF_file = function() {
#   listSNPs = list()
#   for(i in 1:10){
#     mafs = read.csv(sprintf(file.path(homeDirectory, '/data/snp_freqs/chr%02d_maf.csv'), i) , header = T)
#     listSNPs[[i]] = mafs
#   }
#   mafDF = bind_rows(listSNPs)
#   ## choose threshold of maf above 0.05 to prevent inflation
#   mafDF = mafDF %>% filter(x > 0.05)
#   return(mafDF)
# }

### reads Panzea file of maize v2 to v4 mapping
### RETURNS v4_coords
# read_v4_coords_file = function() {
#   v4_coords = read.delim(file.path(homeDirectory, '/data/gbs2.7_agpv4.bed'), header = F)
#   colnames(v4_coords) = c('CHR', 'BP', 'END', 'SNP', 'STRAND')
#   v4_coords$CHR = as.integer(v4_coords$CHR)
#   v4_coords = v4_coords %>% drop_na(CHR)
#   return(v4_coords)
# }


### reads all GEA GWAS p-value results from transformed climate variable data
### RETURNS cumulative results from all chromosomes
read_GEA_results = function(directory, pattern) {
  GEA_results = bind_rows(lapply(list.files(directory, pattern = pattern, full.names = T), 
                                  read.csv))
  GEA_results = GEA_results[, c(3, 1, 2, 4, 9)]
  colnames(GEA_results) = c('SNP', 'CHR', 'BP', 'maf', 'P')
  GEA_results$fdr = p.adjust(GEA_results$P, method = 'BH')
  return(GEA_results)
}

### reads all JOINT GWAS p-value results because formatting is wonky
read_GWAS_results = function(directory, pattern) {
  GWAS_results = bind_rows(lapply(list.files(directory, pattern = pattern, full.names = T), 
                                  read.csv))
  GWAS_results = GWAS_results[, c(1, 2, 3, 6, 10)]
  colnames(GWAS_results) = c('SNP', 'CHR', 'BP', 'maf', 'P')
  GWAS_results$fdr = p.adjust(GWAS_results$P, method = 'BH')
  return(GWAS_results)
}

### reads all GEA GWAS effect size results (beta hats) from transformed climate variable data
### RETURNS cumulative results from all chromosomes
read_GWAS_effects = function(directory, pattern) {
  beta_hats = bind_rows(lapply(list.files(directory, pattern = pattern, full.names = T), 
                               read.csv))
  # beta_hats = cbind(clim_results, beta_hats)
  return(beta_hats)
}

read_GWAS_SEs = function(directory, pattern) {
  SEs = bind_rows(lapply(list.files(directory, pattern = pattern, full.names = T), 
                         read.csv))
  return(SEs)
}

read_GEMMA_results = function(directory, traitN, pattern) {
  gemma_directory = file.path(directory, traitN, 'rep_01')
  gemma_results = bind_rows(lapply(list.files(gemma_directory, pattern = pattern, full.names = T), 
                   read.csv))
  # gemma_results = read.csv(file.path(directory, traitN, 'rep_01', 'p_wald_chr04.txt'))
  gemma_results = separate(data = gemma_results,
                           col = X, into = c('Chr', 'BP'), , sep = '_', remove = F)[c('X', 'Chr', 'BP', traitN)]
  gemma_results$Chr = as.numeric(gsub('S', '', gemma_results$Chr))
  gemma_results$BP = as.numeric(gemma_results$BP)
  return(gemma_results)
}

# ### reads all GEA GWAS p-value results from transformed climate variable data
# ### RETURNS cumulative results from all chromosomes
# read_GEA_results = function() {
#   clim_results = bind_rows(lapply(list.files(file.path(homeDirectory, '/output/clim_log_results'), pattern = 'eGWAS_results', full.names = T), 
#                                   read.csv))
#   clim_results = clim_results[, c(4, 1, 2, 5, 10)]
#   colnames(clim_results) = c('SNP', 'CHR', 'BP', 'maf', 'P')
#   return(clim_results)
# }
# 
# ### reads all GEA GWAS effect size results (beta hats) from transformed climate variable data
# ### RETURNS cumulative results from all chromosomes
# read_GEA_effects = function() {
#   beta_hats = bind_rows(lapply(list.files(file.path(homeDirectory, '/output/clim_log_results'), pattern = 'eGWAS_beta_hats', full.names = T), 
#                                read.csv))
#   # beta_hats = cbind(clim_results, beta_hats)
#   return(beta_hats)
# }

### filters results from a GWAS results dataset based on list of SNPs (intended to be for maf filter)
### RETURNS filtered dataset
filter_SNP_results = function(gwas_results, snpDF) {
  gwas_results = gwas_results %>% filter(SNP %in% rownames(snpDF))
  return(gwas_results)
}

# ### adds v4 coordinates and name to GWAS result dataset
# ### RETURNS same GWAS dataset, with v4 coordinates and name
# add_v4_coords = function(gwas_results, v4_coords) {
#   gwas_results = merge(gwas_results, v4_coords, by = 'SNP')
#   # create v4 SNP name for future analysis with different datasets
#   gwas_results$v4_SNP = paste0('S', gwas_results$CHR, '_', gwas_results$END)
#   return(gwas_results)
# }

## load results from clumping of GEA genes based on LD = 0.3 and specified p-value
## RETURNS dataframe of clumping results
load_clumped = function(pvalue = 1e-05) {
  genome_clumped = read.csv(here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', sprintf('envGWAS_results.genomeclumped_%s.csv', toString(pvalue))))
  # genome_clumped = genome_clumped_03 %>% filter(P < pvalue)
  return(genome_clumped)
}


### function that processes the effect size and p-value results of each phenotype joint GWAS p-value result and effects
### RETURNS combined p-value and effect size dataframe, filtered for maf and calculated p-value
process_GWAS = function(df_results, df_betas) {
  df_results = cbind(df_results, df_betas)
  # df_results = filter_SNP_results(df_results, mafDF)
  df_results$logp = -log10(df_results$P)
  return(df_results)
}

compare_GWAS_effects = function(df, top_hits, gea = F) {
  ### function that takes a Joint result from processJointresults and the top hits from envGWAS
  ### and adds columns that identify top SNP hits as well as samples SNPs with matching maf
  # create factor column from maf using cut (20 intervals in format [0, 0.025))
  if(gea == F) {
    df$maf_factor = cut(df$maf, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
  } else if(gea == T) {
    df$maf_factor = cut(df$maf.gea, breaks = seq(0, 0.50001, 0.025), include.lowest = T, right = F)
  }
  df$is_top_hit = df$SNP %in% top_hits
  
  # use purrr to identify matching maf snps
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
    mutate(non_top_hits = map(data, ~filter(., !is_top_hit))) %>% 
    # create new column by sampling non-top-hits using frequency of top hits
    mutate(samp = map2(non_top_hits, count, sample_n))
  
  # create boolean column for sampled snps that are in same distribution of top hits
  df$matching = df$SNP %in% 
    (nested_df %>% 
       select(samp) %>% 
       unnest(samp) %>% 
       pull(SNP))
  return(df)
}


## merges phenotype joint GWAS with GEA GWAS results
## RETURNS merged dataframe with specific suffixes for maf and P-value
merge_JGWAS_GEA = function(jgwas, gea) {
  # result_df = merge(gea[c('SNP', 'CHR', 'BP', 'maf', 'P', 'altitude..X', 'meanTemp..X', 'annualPrecipitation..X')],
  result_df = merge(gea,
                    jgwas %>% dplyr::select(c(starts_with('Experiment'), 'SNP', 'maf', 'P', 'logp')),
                    by = 'SNP', suffixes = c('.gea', '.jgwas'))
  return(result_df)
}

## polarize JGWAS effect sizes at each trial by the GEA GWAS effect sign, for the passed climate variable.
## i.e. if meanTemp..X is negative at a given SNP, all JGWAS trial effect signs are flipped at that SNP
## thus if that SNP's alternate allele is negative with temperature...
## RETURNS dataframe with polarized values
polarize_JGWAS_GEA = function(jgwas_gea, variable) {
  
  # scale effect sizes by each trial to compare across trials
  jgwas_gea[, grepl('Experimentom', names(jgwas_gea))] = scale(jgwas_gea %>% select(starts_with('Experimentom')))
  
  return(jgwas_gea %>% 
           mutate_at(vars(starts_with('Experimentom')), 
                     funs(ifelse(jgwas_gea[, variable] < 0, -., .))))
}



## polarize and weight JGWAS effect sizes at each trial by multiplying by GEA GWAS effect size, for the passed climate variable.
## thus if that SNP's alternate allele is negative with temperature...
## RETURNS dataframe with extra ".weighted" columns
weight_JGWAS_GEA = function(jgwas_gea, variable) {
  ## this SNP is problematic - see notes
  jgwas_gea = jgwas_gea %>% filter(SNP != 'S9_115491734')
  # scale effect sizes by each trial to compare across trials
  jgwas_gea[, grepl('Experimentom', names(jgwas_gea))] = scale(jgwas_gea %>% select(starts_with('Experimentom')))
  
  GEA_top_hits_sum = jgwas_gea %>% filter(is_top_hit) %>% pull(variable) %>% abs() %>% sum()
  GEA_matching_sum = jgwas_gea %>% filter(matching) %>% pull(variable) %>% abs() %>% sum()
  gwph_weighted = jgwas_gea %>% 
    mutate(
      across(
        .cols = starts_with('Experimentom'),
        .names = '{.col}.weighted',
        ~ .x * jgwas_gea[, variable] / (jgwas_gea$is_top_hit * GEA_top_hits_sum + jgwas_gea$matching * GEA_matching_sum)
      )
    ) 
  
  # %>% select(ends_with('.weighted')) %>% 
  return(gwph_weighted)
}

## calculate average effect sizes at each trial by top_hit/matched and differences between the two
## RETURNS dataframe of differences for each trial
calculateTrialEffectAverages = function(df, trialclim, weighted = F) {
  if(weighted == T) {
    trialEffectAves = df %>% 
      filter(is_top_hit | matching) %>% 
      pivot_longer(cols = ends_with('.weighted'),
                   names_to = 'trial_variables',
                   values_to = 'trial_values') %>% 
      group_by(trial_variables, is_top_hit) %>% 
      summarize(avg = mean(trial_values, na.rm = T))
  } else {
    trialEffectAves = df %>% 
      filter(is_top_hit | matching) %>% 
      pivot_longer(cols = starts_with('Experimentom'),
                   names_to = 'trial_variables',
                   values_to = 'trial_values') %>% 
      group_by(trial_variables, is_top_hit) %>% 
      summarize(avg = mean(trial_values, na.rm = T))
  }
  trialEffectAves$trial_name = str_extract(trialEffectAves$trial_variables, "(?<=Experimento)(.*)(?=.X)")
  trialEffectAves = merge(trialEffectAves, trialclim, by.x = 'trial_name', by.y = 'Experimento')
  trialEffectAves = merge(trialEffectAves, trialinfo[, c('Experimento', 'Estado', 'Trial_elevation')], 
                          by.x = 'trial_name', by.y = 'Experimento')
  return(trialEffectAves)
}

calculateTrialEffectDifferences = function(trialEffectAves) {
  trialEffectDiffs = trialEffectAves %>% 
    group_by(trial_variables) %>% 
    summarise(diff = diff(avg))
  trialEffectDiffs$trial_name = str_extract(trialEffectDiffs$trial_variables, "(?<=Experimento)(.*)(?=.X)")
  trialEffectDiffs = merge(trialEffectDiffs, trialclim, by.x = 'trial_name', by.y = 'Experimento')
  return(trialEffectDiffs)
}


### reads all GxE GWAS p-value results from transformed climate variable data
### RETURNS cumulative results from all chromosomes
read_GxE_GWAS_results = function(directory, pattern) {
  GWAS_results = bind_rows(lapply(list.files(directory, pattern = pattern, full.names = T), 
                                  read_csv))
  GWAS_results = GWAS_results[, c(4, 1, 2, 5, 13)]
  colnames(GWAS_results) = c('SNP', 'CHR', 'BP', 'maf', 'P')
  GWAS_results$fdr = p.adjust(GWAS_results$P, method = 'BH')
  return(GWAS_results)
}