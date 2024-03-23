####################################################################################
# Process results from phenotypic trials
#
# Author: Forrest Li
# Phenotypic prediction analysis
####################################################################################


source(here::here('config.R'))

library(emmeans)
library(multcomp)

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
just_clim = finalMat[-c(1)]
blups_std = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
blups_deregressed = read.csv(here(phenotype_data_dir, 'blups_deregressed.csv'))
envMat = separate(data = finalMat, col = Unique.ID, into = c('SampleID', 'Unique.ID'), sep = ':') %>% subset(select = -c(Unique.ID))

phenotypes_env = blups_std %>% filter(SampleID %in% envMat$SampleID) %>% merge(envMat, by = 'SampleID')

## read phenotypic joint GWAS results and effect sizes for dataset_01_INT which is de-regressed BLUPs
{
  dtf_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_01_INT'), 'GWAS_results')
  dtf_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
  asi_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_01_INT'), 'GWAS_results')
  asi_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
  gwph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_01_INT'), 'GWAS_results')
  gwph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
  ph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_01_INT'), 'GWAS_results')
  ph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
  fw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_01_INT'), 'GWAS_results')
  fw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
  bcw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_01_INT'), 'GWAS_results')
  bcw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_01_INT'), 'GWAS_beta_hats_chr')
}

## done for dataset_03_INT which represents GWAS done on non de-regressed BLUPs (blups_std.csv)
{
  dtf_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_03_INT'), 'GWAS_results')
  dtf_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
  asi_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_03_INT'), 'GWAS_results')
  asi_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
  gwph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_03_INT'), 'GWAS_results')
  gwph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
  ph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_03_INT'), 'GWAS_results')
  ph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
  fw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_03_INT'), 'GWAS_results')
  fw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
  bcw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_03_INT'), 'GWAS_results')
  bcw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
}

## combine joint GWAS results to effect sizes, get -log10 on p-values
{
  dtf_process = process_GWAS(dtf_results, dtf_betas)
  bcw_process = process_GWAS(bcw_results, bcw_betas)
  asi_process = process_GWAS(asi_results, asi_betas)
  gwph_process = process_GWAS(gwph_results, gwph_betas)
  fw_process = process_GWAS(fw_results, fw_betas)
  ph_process = process_GWAS(ph_results, ph_betas)
}

## import SE data
{
  dtf_SE = process_GWAS(dtf_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_03_INT'), 'GWAS_SEs'))
  asi_SE = process_GWAS(asi_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_03_INT'), 'GWAS_SEs'))
  gwph_SE = process_GWAS(gwph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_03_INT'), 'GWAS_SEs'))
  ph_SE = process_GWAS(ph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_03_INT'), 'GWAS_SEs'))
  fw_SE = process_GWAS(fw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_03_INT'), 'GWAS_SEs'))
  bcw_SE = process_GWAS(bcw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_03_INT'), 'GWAS_SEs'))
  
}

####################################################################################
## write a table for seeing how many accessions are in each trial for each trait
####################################################################################
phenotypes = c('GrainWeightPerHectareCorrected', 'PlantHeight', 'ASI', 'FieldWeight', 'BareCobWeight', 'DaysToFlowering')
trial_trait_table = setNames(data.frame(matrix(ncol = 6, nrow = 31), 
                                        row.names = blups_std$Experimento %>% unique()),
                             phenotypes)
for(trait in phenotypes) {
  for(experiment in blups_std$Experimento %>% unique()) {
    num_acc = phenotypes_env %>% filter(Trait == trait, Experimento == experiment) %>% pull(SampleID) %>% unique() %>% length()
    trial_trait_table[experiment, trait] = num_acc
  }
}

write.csv(trial_trait_table, here(phenotype_data_dir, 'accessions_per_trial_per_trait.csv'))

## do this for deregressed trials, as some are dropped
trial_trait_table_deregressed = setNames(data.frame(matrix(ncol = 6, nrow = 31), 
                                        row.names = blups_deregressed$Experimento %>% unique()),
                             phenotypes)
for(trait in phenotypes) {
  for(experiment in blups_deregressed$Experimento %>% unique()) {
    num_acc = phenotypes_env %>% filter(Trait == trait, Experimento == experiment) %>% pull(SampleID) %>% unique() %>% length()
    trial_trait_table_deregressed[experiment, trait] = num_acc
  }
}

write.csv(trial_trait_table_deregressed, here(phenotype_data_dir, 'accessions_per_trial_per_trait_deregressed.csv'))



####################################################################################
## Just print out samples of matched SNPs and be done with this farce
####################################################################################
gea_matching_samples = as.data.frame(replicate(20, compare_GWAS_effects(gea_results, top_hits_clumped, gea = T) %>% filter(matching == T) %>% pull(SNP)))
write.table(gea_matching_samples, here(genetic_data_dir, 'gea_matching_sampled_SNPs.txt'), col.names = F, row.names = F)

## print out GEA SNPs under 1e-2 p-value for plink clumping and LD analysis
SNP_LD = gea_results %>% filter(P < 1e-2) %>% select(SNP)
write.csv(SNP_LD, 'Analyses/GEA_output/multivariate_results/GEA_SNP_1e2_forplink.txt', row.names = F, quote = F)

plink_analysis = fread('Analyses/GEA_output/multivariate_results/clumped/plink/LD-analysis.ld')

manual_manhattan(gea_results, plink_analysis$SNP_B, chr_filter = 7)

####################################################################################
var_lookup = c(tmin = 'tmin..X', tmax = 'tmax..X', trange = 'trange..X', precipTot = 'precipTot..X', 
               aridityMean = 'aridityMean..X', rhMean = 'rhMean..X', elevation = 'elevation..X')
gea_fx_long = gea_beta_hats %>% 
  rename(all_of(var_lookup)) %>% 
  pivot_longer(cols = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation'), 
               names_to = 'env_var')

all_traits_fx = do.call(rbind, purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process), 
                                           list('bare cob weight', 'grain weight per hectare', 'field weight', 'plant height', 'days to flowering', 'ASI'),
                                           get_heat_map))
####################################################################################
## phenotypic prediction
####################################################################################

trial_info = read.csv('Phenotype_data/Trial_info.csv')
blups_std = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
blups_deregressed = read.csv(here(phenotype_data_dir, 'blups_deregressed.csv'))

accessions_per_trial = read.csv(here(phenotype_data_dir, 'accessions_per_trial_per_trait.csv'))
trial_worldclim = read.csv('Phenotype_data/Trial_worldclim.csv')
trial_master = merge(trial_worldclim, trial_info, by = 'Experimento')

## plot correlations between trials for all raw worldclim variables
cor_between_trials = data.frame(cor(t(trial_worldclim[4:22])), 
                                trial_worldclim$Experimento)
colnames(cor_between_trials) = c(trial_worldclim$Experimento, 'Trial')

## make tables for prediction ability values
all_traits_prediction_ability_deregressed = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
  cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups/modelPrediction_%s_results_resid.csv', trait_i)), trait_i)
})))

all_traits_prediction_ability_std = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
  cbind(fread(sprintf('Analyses/PhenotypicPrediction/fiveModels/modelPrediction_%s_results_resid.csv', trait_i)), trait_i)
})))

write.csv(all_traits_prediction_ability_deregressed, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_deregressed.csv'))
write.csv(all_traits_prediction_ability_std, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_std.csv'))

## analyze mean prediction ability across trials
all_traits_prediction_ability_deregressed_mean = all_traits_prediction_ability_deregressed %>% group_by(trait_i) %>% summarise_at(vars(c('rf_env', 'lm_PC5', 'lm_all_SNPs', 
                                                                                                     'lm_enriched_SNPs', 'lm_enriched_SNPs', 'lm_matching_SNPs')), 
                                                                                                   .funs = c(mean = 'mean'))
write.csv(all_traits_prediction_ability_deregressed_mean, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_deregressed_mean.csv'))


all_traits_prediction_ability_deregressed_melt = pivot_longer(all_traits_prediction_ability_deregressed, 
                                                              cols = c('rf_env', 'lm_PC5', 'lm_all_SNPs', 'lm_enriched_SNPs', 'lm_enriched_SNPs', 'lm_matching_SNPs'),
                                                              names_to = 'model',
                                                              values_to = 'r')
## across all traits
{
  prediction_comparison_model = lm(r ~ Experimento + fold + trait_i + model, 
                                   data = all_traits_prediction_ability_deregressed_melt)
  anova(prediction_comparison_model)
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1 - 0.05/5)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/5, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value*5)
  prediction_comparison_effects_summary
}

## just for field weight
{
  summary(all_traits_prediction_ability_deregressed %>% filter(trait_i == 'FieldWeight'))
  prediction_comparison_model = lm(r ~ Experimento + fold + model, 
                                   data = all_traits_prediction_ability_deregressed_melt %>% filter(trait_i == 'FieldWeight'))
  anova(prediction_comparison_model)
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1 - 0.05/1)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/1, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value*1)
  prediction_comparison_effects_summary
}