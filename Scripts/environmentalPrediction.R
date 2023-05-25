####################################################################################
# Process MegaLMM environmental prediction results
#
# Author: Forrest Li
# Env Prediction analysis
####################################################################################

source(here::here('config.R'))

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

# correlation plot for just withheld data
ggplot(observed_comparison %>% filter(isWithheld & climate_variables != 'rhMean'), aes(climate_values.x, climate_values.y, color = factor(isWithheld))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = c('red')) +
  facet_wrap('climate_variables', scales = 'free', labeller = labeller(climate_variables = c())) +
  xlab('Observed climate variable ranking') +
  ylab('Predicted climate variable ranking') +
  # labs(color = 'Predicted accession') +
  stat_cor(aes(label = ..rr.label..), color = 'blue', geom = 'label') +
  # ggtitle('Prediction correlation') +
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20)) 
