####################################################################################
# Process MegaLMM environmental prediction results
#
# Author: Forrest Li
# Env Prediction analysis
####################################################################################

source(here::here('config.R'))
library(doParallel)
library(foreach)

G_clim_pred = as.matrix(read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_01', 'G_hat-rep_01.csv'), row.names = 1))
R_clim_pred = as.matrix(read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_01', 'R_hat-rep_01.csv'), row.names = 1))
G_clim = as.matrix(read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-rep_02', 'G_hat.csv'), row.names = 1))
R_clim = as.matrix(read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-rep_02', 'R_hat.csv'), row.names = 1))

U_hat = read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_01', 'U_hat-rep_01.csv'), row.names = 1)
withheld = read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_01', 'withheld-accessions-rep_01.csv'), header = T)

U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
U_hat = U_hat[match(finalMat$Unique.ID, U_hat$SampleID),]
U_hat$isWithheld = U_hat$SampleID %in% withheld$validationIDs

variables = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
colnames(finalMat) = c('SampleID', variables)

observed_comparison = merge(U_hat %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'climate_values'
                              ),
                            finalMat %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'climate_values'
                              ), by = c('SampleID', 'climate_variables'))

ggplot(observed_comparison, aes(climate_values.x, climate_values.y, color = factor(isWithheld))) +
  geom_point() +
  scale_color_manual(labels = c('training', 'predicted'), values = c('black', 'red')) +
  facet_wrap('climate_variables', scales = 'free') +
  xlab('Observed climate variables') +
  ylab('Predicted U_hat') +
  labs(color = 'Predicted accession') +
  # stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
  ggtitle("Prediction correlation") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5), text = element_text(size = 14))

# correlation plot for just withheld data - make sure to swap x-y axes
(environmental_prediction = ggplot(observed_comparison %>% filter(isWithheld & climate_variables != 'rhMean'), aes(climate_values.y, climate_values.x, color = factor(isWithheld))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = c('red')) +
  facet_wrap('climate_variables', scales = 'free', labeller = labeller(climate_variables = c())) +
  xlab('Predicted climate variable ranking') +
  ylab('Observed climate variable ranking') +
  # labs(color = 'Predicted accession') +
  stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
  # ggtitle('Prediction correlation') +
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20)) 
)

###################################################################################################
## perform for k-fold testing
###################################################################################################
results = list()
for(i in 1:5) {
  U_hat = read.csv(here::here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', sprintf('U_hat-rep_05-fold_%02d.csv', i)), row.names = 1)
  withheld = read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', sprintf('withheld-accessions-rep_05-fold_%02d.csv', i)), header = T)
  U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)
  U_hat = U_hat %>% filter(SampleID %in% withheld$validationIDs)
  results[[i]] = U_hat
}

combined = bind_rows(results)

# combined = read.csv(here::here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', 'U_hat-rep_05-combined.csv'), row.names = 1)
# U_hat = read.csv(here::here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', 'U_hat-rep_05-fold_05.csv'), row.names = 1)
# withheld = read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', 'withheld-accessions-rep_05-fold_05.csv'), header = T)
# 
# U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)
# 
finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
variables = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')
U_hat = U_hat[match(finalMat$Unique.ID, U_hat$SampleID),]
colnames(finalMat) = c('SampleID', variables)
GEA_enriched = read.csv(here('Genetic_data/Imputed_V4/GEA_clumped_SNPs_genotype.csv'))
finalMat = merge(finalMat, GEA_enriched, by.x = 'SampleID', by.y = 'X')


observed_comparison = merge(combined %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'climate_values'
                              ),
                            finalMat %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'climate_values'
                              ), by = c('SampleID', 'climate_variables'))


# correlation plot for just withheld data - make sure to swap x-y axes
(environmental_prediction = ggplot(observed_comparison %>% filter(climate_variables != 'rhMean'), 
                                   aes(climate_values.x, climate_values.y)) +
    geom_point(show.legend = F) +
    facet_wrap('climate_variables', scales = 'free', labeller = labeller(climate_variables = c())) +
    xlab('Predicted climate variable ranking') +
    ylab('Observed climate variable ranking') +
    # labs(color = 'Predicted accession') +
    stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
    # ggtitle('Prediction correlation') +
    theme(strip.text.x = element_text(size = 20),
          axis.title = element_text(size = 20)) 
)

(enviromental_prediction_allele = ggplot(observed_comparison %>% filter(climate_variables != 'rhMean'), 
                                   aes(climate_values.x, climate_values.y,
                                       color = S9_148365695)) +
    geom_point(size = 0.5) +
    scale_color_continuous(low = 'red', high = 'blue') +
    facet_wrap('climate_variables', scales = 'free', labeller = labeller(climate_variables = c())) +
    xlab('Predicted climate variable ranking') +
    ylab('Observed climate variable ranking') +
    # labs(color = 'Predicted accession') +
    stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
    # ggtitle('Prediction correlation') +
    theme(strip.text.x = element_text(size = 20),
          axis.title = element_text(size = 20)) 
)


###################################################################################################
## perform for spatial fold testing
###################################################################################################
spatial_results = list()
rep = 2
for(i in 1:10) {
  U_hat = read.csv(here::here(analyses_dir, 'MegaLMM_output', sprintf('dataset_clim-spatial_prediction-rep_%02d/U_hat-rep_%02d-fold_%02d.csv', rep, rep, i)), row.names = 1)
  withheld = read.csv(here(analyses_dir, 'MegaLMM_output', sprintf('dataset_clim-spatial_prediction-rep_%02d/withheld-accessions-rep_%02d-fold_%02d.csv', rep, rep, i)), header = T)
  U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)
  U_hat = U_hat %>% filter(SampleID %in% withheld$validationIDs)
  spatial_results[[i]] = U_hat
}

spatial_combined = bind_rows(spatial_results, .id = 'fold')

spatial_comparison = merge(spatial_combined %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'predicted_values'
                              ),
                            finalMat %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'observed_values'
                              ), by = c('SampleID', 'climate_variables'))

(environmental_prediction_spatial = ggplot(spatial_comparison %>% filter(climate_variables != 'rhMean'), 
                                   aes(predicted_values, observed_values, color = fold)) +
    geom_point(size = 0.3) +
    facet_wrap('climate_variables', scales = 'free', labeller = labeller(climate_variables = c())) +
    xlab('Predicted climate variable ranking') +
    ylab('Observed climate variable ranking') +
    # labs(color = 'Predicted accession') +
    stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
    # ggtitle('Prediction correlation') +
    theme(strip.text.x = element_text(size = 20),
          axis.title = element_text(size = 20)) 
)

random_results = list()
rep = 3
for(i in 1:10) {
  U_hat = read.csv(here::here(analyses_dir, 'MegaLMM_output', sprintf('dataset_clim-spatial_prediction-rep_%02d/U_hat-rep_%02d-fold_%02d.csv', rep, rep, i)), row.names = 1)
  withheld = read.csv(here(analyses_dir, 'MegaLMM_output', sprintf('dataset_clim-spatial_prediction-rep_%02d/withheld-accessions-rep_%02d-fold_%02d.csv', rep, rep, i)), header = T)
  U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)
  U_hat = U_hat %>% filter(SampleID %in% withheld$validationIDs)
  random_results[[i]] = U_hat
}

random_combined = bind_rows(random_results, .id = 'fold')

all_combined = bind_rows(spatial_combined %>% 
                           mutate(fold_type = 'spatial') %>% 
                             pivot_longer(
                               cols = variables,
                               names_to = 'climate_variables',
                               values_to = 'predicted_values'), 
                           random_combined %>% 
                           mutate(fold_type = 'random') %>% 
                           pivot_longer(
                                    cols = variables,
                                    names_to = 'climate_variables',
                                    values_to = 'predicted_values')
) %>% merge(finalMat %>% 
              pivot_longer(
                cols = variables,
                names_to = 'climate_variables',
                values_to = 'observed_values'), 
            by = c('SampleID', 'climate_variables'))

res_spclust = foreach(trait=variables, .combine = rbind) %:% 
  foreach(i_fold_type = c('spatial', 'random'), .combine = rbind) %:%
    foreach(i_fold = 1:10, .combine = rbind) %do% {
      observed = all_combined %>% filter(climate_variables == trait, fold == i_fold, fold_type == i_fold_type) %>% pull(observed_values)
      predicted = all_combined %>% filter(climate_variables == trait, fold == i_fold, fold_type == i_fold_type) %>% pull(predicted_values)
      data.frame(fold = i_fold, trait = trait, fold_type = i_fold_type, r = cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2)))
  }


## plots for correlation and RMSE predictive accuracy for both spatial and random sampling
(r_boxplot = ggplot(res_spclust,aes(x=trait,y=r, color = fold_type)) + 
  geom_boxplot() +
    expand_limits(y = c(0, 1)) +
  ggtitle('Pearson r correlation for 10 folds')
)
(rmse_boxplot = ggplot(res_spclust,aes(x=trait,y=rmse, color = fold_type)) + 
  geom_boxplot() + 
  expand_limits(y = c(0, 1)) +
  ggtitle('RMSE predictive accuracy for 10 folds')
)

## plots for correlation and RMSE spatial and random per fold for different traits
(spatial_R = ggplot(subset(res_spclust,fold_type == 'spatial'),aes(x=trait,y=r)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  ggtitle('Pearson r correlation for spatial sampling') +
  expand_limits(y=c(0,1))
)

(spatial_RMSE = ggplot(subset(res_spclust,fold_type == 'spatial'),aes(x=trait,y=rmse)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  ggtitle('RMSE predictive accuracy for spatial sampling') +
  expand_limits(y=c(0,1))
)

(random_R = ggplot(subset(res_spclust,fold_type == 'random'),aes(x=trait,y=r)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  ggtitle('Pearson r correlation for random sampling') +
    expand_limits(y=c(0,1))
)

(random_RMSE = ggplot(subset(res_spclust,fold_type == 'random'),aes(x=trait,y=rmse)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  ggtitle('RMSE predictive accuracy for random sampling') +
    expand_limits(y=c(0,1))
)

(figure4 = ggarrange(
  environmental_prediction,
  ggarrange(r_boxplot, rmse_boxplot, ncol = 2, labels = c('A', 'B'), common.legend = T, legend = 'right'), 
  ggarrange(spatial_R, random_R, spatial_RMSE, random_RMSE, ncol = 2, nrow = 2, labels = c('C', 'D', 'E', 'F'), common.legend = TRUE, legend = 'right'), nrow = 3))
png(here(plot_dir, 'Manuscript', 'prediction-of-environment.png'), width = 982, height = 982)
print(figure4)
dev.off()    
    



