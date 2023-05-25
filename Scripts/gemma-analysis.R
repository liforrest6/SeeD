####################################################################################
# Process envGWAS results from univariate GEA analyses
#
# Author: Forrest Li
# EnvGWAS analysis for GEMMA results
####################################################################################

source(here::here('config.R'))

clim_results = read_GEA_results(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_results')

gemma_dir = here('Analyses', 'GEA_output', 'GEMMA_univariate')
traitN = 'precipTot'
gemma_results = read.csv(file.path(gemma_dir, 'precipTot', 'rep_01', 'p_wald_chr04.txt'))
gemma_results = separate(data = gemma_results, 
                         col = X, into = c('Chr', 'BP'), , sep = '_', remove = F)[c('X', 'Chr', 'BP', traitN)]
gemma_results$Chr = as.numeric(gsub('S', '', gemma_results$Chr))
gemma_results$BP = as.numeric(gemma_results$BP)

{
  precipTot_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'precipTot', 'p_wald')
  tmin_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmin', 'p_wald')
  tmax_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmax', 'p_wald')
  trange_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'trange', 'p_wald')
  rhMean_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'rhMean', 'p_wald')
  aridityMean_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'aridityMean', 'p_wald')
  elevation_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'elevation', 'p_wald')
}

