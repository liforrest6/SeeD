####################################################################################
# Process envGWAS results from univariate GEA analyses
#
# Author: Forrest Li
# EnvGWAS analysis for GEMMA results
####################################################################################

source(here::here('config.R'))

library('ggVennDiagram')

gea_results = read_GEA_results(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_results')

gemma_dir = here('Analyses', 'GEA_output', 'GEMMA_univariate')
# traitN = 'precipTot'
# gemma_results = read.csv(file.path(gemma_dir, 'precipTot', 'rep_01', 'p_wald_chr04.txt'))
# gemma_results = separate(data = gemma_results, 
#                          col = X, into = c('Chr', 'BP'), , sep = '_', remove = F)[c('X', 'Chr', 'BP', traitN)]
# gemma_results$Chr = as.numeric(gsub('S', '', gemma_results$Chr))
# gemma_results$BP = as.numeric(gemma_results$BP)

{
  precipTot_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'precipTot', 'p_wald')
  tmin_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmin', 'p_wald')
  tmax_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmax', 'p_wald')
  trange_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'trange', 'p_wald')
  rhMean_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'rhMean', 'p_wald')
  aridityMean_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'aridityMean', 'p_wald')
  elevation_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'elevation', 'p_wald')
}

gea_results_top = gea_results[gea_results$P < 1e-6, 'SNP']
tmax_results_top = tmax_results[tmax_results$tmax < 1e-5, 'V1']
precipTot_results_top = precipTot_results[precipTot_results$precipTot < 1e-5, 'V1']
elevation_results_top = elevation_results[elevation_results$elevation < 1e-5, 'V1']
trange_results_top = trange_results[trange_results$trange < 1e-5, 'V1']
tmin_results_top = tmin_results[tmin_results$tmin < 1e-5, 'V1']
rhMean_results_top = rhMean_results[rhMean_results$rhMean < 1e-5, 'V1']
aridityMean_results_top = aridityMean_results[aridityMean_results$aridityMean < 1e-5, 'V1']


all_gwas_top_hits = list(GEA = gea_results_top,
                         tmax = tmax_results_top,
                         precipTot = precipTot_results_top,
                         elevation = elevation_results_top,
                         aridityMean = aridityMean_results_top,
                         tmin = tmin_results_top,
                         rhMean = rhMean_results_top,
                         trange = trange_results_top)

ggVennDiagram(all_gwas_top_hits)
