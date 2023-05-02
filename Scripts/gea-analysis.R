####################################################################################
# Process envGWAS results from multivariate GEA analyses
#
# Author: Forrest Li
# EnvGWAS analysis
####################################################################################

source(here::here('config.R'))

clim_results = read_GWAS_results(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_results')

top_hits_clumped = load_clumped()

qqman::qq(clim_results$P, main = 'qqplot, maf > 0.01')

qqman::qq(clim_results %>% filter(maf > 0.05) %>% pull(P), main = 'qqplot, maf > 0.05')

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
just_clim = finalMat[-c(1)]
library(randomForest)
forest = randomForest(elevation ~ ., data = just_clim, ntree=1000,
             keep.forest=FALSE, importance=TRUE)

