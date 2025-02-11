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

##### perform for random k-fold testing################################################################################

environmental_prediction_results = list()
for(i in 1:5) {
  U_hat = read.csv(here::here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', sprintf('U_hat-rep_05-fold_%02d.csv', i)), row.names = 1)
  withheld = read.csv(here(analyses_dir, 'MegaLMM_output', 'dataset_clim-prediction-rep_05', sprintf('withheld-accessions-rep_05-fold_%02d.csv', i)), header = T)
  U_hat$SampleID = str_split_i(rownames(U_hat), '::', 1)
  U_hat = U_hat %>% filter(SampleID %in% withheld$validationIDs)
  environmental_prediction_results[[i]] = U_hat
}

combined = bind_rows(environmental_prediction_results)

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
finalMat_genotype = merge(finalMat, GEA_enriched, by.x = 'SampleID', by.y = 'X')


observed_comparison = merge(combined %>% 
                              pivot_longer(
                                cols = variables,
                                names_to = 'climate_variables',
                                values_to = 'climate_values'
                              ),
                            finalMat_genotype %>% 
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
    geom_point(size = 0.3) +
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


##### perform for spatial fold testing  ###################################################################

## read all spatial folds
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

## merge with real climate values
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

## plot spatial accuracy points
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

## read all random folds
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

## combine both spatial folds and random folds into one dataframe
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



##### calculate means for each trait and fold_type ############################################################################
## check if there is sig difference in means between fold types

## generate dataframe looking at aggregate accuracies
res_spclust = foreach(trait=variables, .combine = rbind) %:% 
  foreach(i_fold_type = c('spatial', 'random'), .combine = rbind) %:%
  foreach(i_fold = 1:10, .combine = rbind) %do% {
    observed = all_combined %>% filter(climate_variables == trait, fold == i_fold, fold_type == i_fold_type) %>% pull(observed_values)
    predicted = all_combined %>% filter(climate_variables == trait, fold == i_fold, fold_type == i_fold_type) %>% pull(predicted_values)
    data.frame(fold = i_fold, trait = trait, fold_type = i_fold_type, r = cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2)))
  }

## run ANOVA to see if fold_type affects predictive accuracy
spatial_vs_random_model = lm(r ~ trait + trait:fold_type + fold_type, res_spclust)
spatial_vs_random_emmeans = emmeans(spatial_vs_random_model, by = 'trait', specs = 'fold_type')
anova(spatial_vs_random_model)
spatial_vs_random_emmeans

## calculate aggregate means
res_spclust %>% 
  group_by(trait, fold_type) %>% 
  summarise(mean_r = mean(r), min_r = min(r), max_r = max(r),
            mean_rmse = mean(rmse))

##### make plots for GPoE analyses #####################################################################


## boxplots for correlation and RMSE predictive accuracy for both spatial and random sampling
(r_boxplot = ggplot(res_spclust,aes(x=trait,y=r, color = fold_type)) + 
  geom_boxplot() +
    expand_limits(y = c(0, 1)) +
    scale_color_discrete(name = "Fold type") +
    ylab("Pearson's correlation (r)") +
    xlab('') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -45, size = 5, vjust = rotation_vjust(-45), hjust = rotation_hjust(45)))+
    # labs(color = '10-fold cross-validation')+

    # theme(text = element_text(size = 16),
          # legend.position = 'top') +
    # theme_bw() +
    # theme(axis.text.x = element_text(size = 16),
    #       axis.text.y = element_text(size = 16)) +
  # ggtitle("Pearson's correlation (r) predictive ability")
    ggtitle("")
)

(rmse_boxplot = ggplot(res_spclust,aes(x=trait,y=rmse, color = fold_type)) + 
  geom_boxplot() + 
    theme_bw() +
    scale_color_discrete(name = "Fold type") +
    ylab("RMSE") +
    xlab('') +
    theme(axis.text.x = element_text(angle = -45, size = 5, vjust = rotation_vjust(-45), hjust = rotation_hjust(45))) +
  expand_limits(y = c(0, 1)) +
  # ggtitle('RMSE predictive ability')
    ggtitle('')
)

# ## random CV for GPoE for slides
# (r_boxplot_figure = ggplot(res_spclust ,aes(x=trait,y=r, color = fold_type)) +
#     geom_boxplot() +
#     expand_limits(y = c(0, 1)) +
#     ylab("Pearson's correlation (r)") +
#     xlab('Environmental covariate') +
#     labs(color = '10-fold cross-validation')+
#     theme_bw() +
#     theme(text = element_text(size = 16),
#           legend.position = 'top') +
#     # theme_bw() +
#     # theme(axis.text.x = element_text(size = 16),
#     #       axis.text.y = element_text(size = 16)) +
#     scale_color_manual(values = c('#8ACE00', 'purple'))
#     # ggtitle('Genomic Prediction of Environment (GPoE) predictive ability')
# )


## plots for correlation and RMSE spatial and random per fold for different traits
(spatial_R = ggplot(subset(res_spclust,fold_type == 'spatial'),aes(x=trait,y=r)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  
    # ggtitle("Spatial fold Pearson's r") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -70, size = 8, vjust = rotation_vjust(-70), hjust = rotation_hjust(70)),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) + 
    xlab('') +
    ylab("Pearson's r") +
    labs(color = 'Fold #')+
  expand_limits(y=c(-.5,1))
)

(spatial_RMSE = ggplot(subset(res_spclust,fold_type == 'spatial'),aes(x=trait,y=rmse)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) +  
    # ggtitle("Spatial fold RMSE") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = -70, size = 8, vjust = rotation_vjust(-70), hjust = rotation_hjust(70)),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) + 
    xlab('') +
    ylab("RMSE") +
    labs(color = 'Fold #')+
  expand_limits(y=c(0,1.5))
)

(random_R = ggplot(subset(res_spclust,fold_type == 'random'),aes(x=trait,y=r)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) + 
    # ggtitle("Random fold Pearson's r") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = -70, size = 8, vjust = rotation_vjust(-70), hjust = rotation_hjust(70)),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) + 
    labs(color = 'Fold #')+
    xlab('') +
    ylab("Pearson's r") +
    expand_limits(y=c(-.5,1))
)

(random_RMSE = ggplot(subset(res_spclust,fold_type == 'random'),aes(x=trait,y=rmse)) +
  geom_line(aes(group = interaction(fold,fold_type),color=as.factor(fold))) +
  geom_text(aes(label=fold,group = interaction(fold,fold_type),color=as.factor(fold)),
            size=3) + 
    # ggtitle("Random fold RMSE") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = -70, size = 8, vjust = rotation_vjust(-70), hjust = rotation_hjust(70)),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) + 
    labs(color = 'Fold #')+
    xlab('') +
    ylab("RMSE") +
    expand_limits(y=c(0,1.5))
)

## create map of spatial folds
# read non-transformed data for lat-long coordinates
finalMat_nontransform = read.csv('Env_data/GEA-climate-nontransformed.csv') %>% dplyr::select(c('Unique.ID', 'LatNew', 'LongNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'))
spatial_combined_coords = merge(spatial_combined[c('fold', 'SampleID')], finalMat_nontransform, by.x = 'SampleID', by.y = 'Unique.ID')
spatial_combined_coords$fold = factor(spatial_combined_coords$fold, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
(spatialsample_plot = ggplot(spatial_combined_coords, aes(x = LongNew, y = LatNew, color = fold)) +
    geom_text(aes(label = fold), size = 2) +
    theme_bw() +
    guides(color = 'none') +
    # labs(color = 'Fold #')+
    xlab('Longitude') +
    ylab('Latitude')
)

## manuscript Figure 2
(prediction_of_environment = ggarrange(r_boxplot, rmse_boxplot, ncol = 2, 
                                       labels = c('A', 'B'), common.legend = T, legend = 'bottom'))
(spatial_sampling = ggarrange(spatialsample_plot,
                              ggarrange(spatial_R + theme(axis.text.x = element_blank()), 
                                        random_R + theme(axis.text.x = element_blank()), 
                                        spatial_RMSE, 
                                        random_RMSE, ncol = 2, nrow = 2, labels = c('B', 'C', 'D', 'E'), 
                                        common.legend = TRUE, legend = 'right', heights = c(2, 3)), 
                              labels = c('A'), ncol = 2, widths = c(2, 3)))
  

ggsave('prediction-of-environment.tif',
       plot = prediction_of_environment,
       device = 'tiff',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 80,
       units = 'mm'
)

ggsave('spatial-sampling.pdf',
       plot = spatial_sampling,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 80,
       units = 'mm'
)


spatial_melt_values = spatial_combined_coords %>% pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
                                     names_to = 'env_trait', values_to = 'env_trait_value')
ggplot(spatial_melt_values, aes(y = env_trait_value, x = fold)) +
  geom_violin() +
  facet_wrap('env_trait',scales = 'free')


# (validation_fold_plot = autoplot(clusters_spatial[[1]]))
# png(here(plot_dir, 'Manuscript', 'validation_buffer.png'), width = 400, height = 700)
# print(validation_fold_plot)
# dev.off()    
##### test plot RMSE of certain folds ########################


spatial_coords = spatial_combined[c('SampleID', 'fold')]

acc_ids = geno_sf$V1
clusters_spatial$splits[[1]] %>% analysis()
folds = data.frame(fold =  1:10)
folds$tmin = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(tmin) %>% sd()}) %>% unlist()
folds$tmax = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(tmax) %>% sd()}) %>% unlist()
folds$trange = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(trange) %>% sd()}) %>% unlist()
folds$precipTot = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(precipTot) %>% sd()}) %>% unlist()
folds$aridityMean = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(aridityMean) %>% sd()}) %>% unlist()
folds$rhMean = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(rhMean) %>% sd()}) %>% unlist()
folds$elevation = lapply(c(1:10), function(x){finalMat %>% 
    filter(SampleID %in% (clusters_spatial$splits[[x]] %>% assessment() %>% pull(V1))) %>% 
    pull(elevation) %>% sd()}) %>% unlist()
folds_pivot = folds %>% pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
                                     names_to = 'trait', values_to = 'trait_sd')

folds_sd_rmse = merge(folds_pivot, res_spclust %>% filter(fold_type == 'spatial'), by = c('fold', 'trait'))
ggplot(folds_sd_rmse, aes(x = trait_sd, y = rmse, color = as.factor(fold), label = as.factor(fold))) +
  # geom_point() +
  geom_text(aes(label=as.factor(fold))) +
  facet_wrap('trait')

folds10 = clusters_spatial$splits[[10]] %>% assessment() %>% pull(V1)
folds2 = clusters_spatial$splits[[2]] %>% assessment() %>% pull(V1)

finalMat %>% filter(SampleID %in% folds10) %>% dplyr::select(c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')) %>% apply(., 2, sd)
  
finalMat_nontransform
finalMat_nontransform %>% filter(Unique.ID %in% folds10) %>% dplyr::select(c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')) %>% apply(., 2, sd)
finalMat_nontransform %>% filter(Unique.ID %in% folds2) %>% dplyr::select(c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')) %>% apply(., 2, sd)


lapply(c(1:10), function(x){clusters_spatial$splits[[x]] %>% assessment() %>% as.data.frame() %>% dplyr::select(c('tmin.x', 'tmax.x', 'trange.x', 'precipTot.x', 'aridityMean.x', 'rhMean.x', 'elevation.x')) %>% apply(., 2, sd)})


