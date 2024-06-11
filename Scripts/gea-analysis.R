####################################################################################
# Process results from multivariate GEA analyses
#
# Author: Forrest Li
# GEA analysis
####################################################################################

source(here::here('config.R'))
library(randomForest)

gea_results = read_GEA_results(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_results')
gea_beta_hats = read_GWAS_effects(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_beta_hats')
gea_beta_hats = cbind(gea_results, gea_beta_hats)

clumped_analysis = load_clumped()
top_hits_clumped = load_clumped()$SNP


