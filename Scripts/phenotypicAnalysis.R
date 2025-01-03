####################################################################################
# Process results from phenotypic trials
#
# Author: Forrest Li
# Phenotypic prediction analysis
####################################################################################


source(here::here('config.R'))

library(emmeans)
library(multcomp)
library(lmerTest)

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
just_clim = finalMat[-c(1)]
blups_std = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
blups_deregressed = read.csv(here(phenotype_data_dir, 'blups_deregressed.csv'))
trial_info = read.csv('Phenotype_data/Trial_info.csv')

finalMat$SampleID = separate(data = finalMat, col = Unique.ID, into = c('SampleID', 'Unique.ID'), sep = ':')$SampleID
# envMat = separate(data = finalMat, col = SampleID, into = c('SampleID', 'Unique.ID'), sep = ':') %>% subset(select = -c(Unique.ID))

# phenotypes_env = blups_deregressed %>% filter(SampleID %in% envMat$SampleID) %>% merge(envMat, by = 'SampleID')

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
# {
#   dtf_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_03_INT'), 'GWAS_results')
#   dtf_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
#   asi_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_03_INT'), 'GWAS_results')
#   asi_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
#   gwph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_03_INT'), 'GWAS_results')
#   gwph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
#   ph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_03_INT'), 'GWAS_results')
#   ph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
#   fw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_03_INT'), 'GWAS_results')
#   fw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
#   bcw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_03_INT'), 'GWAS_results')
#   bcw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_03_INT'), 'GWAS_beta_hats_chr')
# }

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
  dtf_SE = process_GWAS(dtf_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'dataset_01_INT'), 'GWAS_SEs'))
  asi_SE = process_GWAS(asi_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'ASI', 'dataset_01_INT'), 'GWAS_SEs'))
  gwph_SE = process_GWAS(gwph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'dataset_01_INT'), 'GWAS_SEs'))
  ph_SE = process_GWAS(ph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'dataset_01_INT'), 'GWAS_SEs'))
  fw_SE = process_GWAS(fw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'dataset_01_INT'), 'GWAS_SEs'))
  bcw_SE = process_GWAS(bcw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'dataset_01_INT'), 'GWAS_SEs'))
  
}
####################################################################################
## add passport race data to prediction results
####################################################################################
passport_data = read.csv('Phenotype_data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt', sep = '\t')
finalMat_race = finalMat %>% merge(passport_data, by.x = 'SampleID', by.y = 'Sample.ID.of.DNA.from.single.plants.used.in.GWAS', all.x = TRUE) %>% dplyr::select(Unique.ID, PrimaryRace)
finalMat_race$PrimaryRace[finalMat_race$PrimaryRace == ''] = 'UNKNOWN'
finalMat_race$PrimaryRace[is.na(finalMat_race$PrimaryRace)] = 'UNKNOWN'
finalMat_race = finalMat_race[match(finalMat$Unique.ID, finalMat_race$Unique.ID),]
# write.csv(finalMat_race, here(env_data_dir, 'GEA-race-classifications.csv'), row.names = F)
 
####################################################################################
## write a table for seeing how many accessions are in each trial for each trait
####################################################################################
traits = c('ASI', 'GrainWeightPerHectareCorrected', 'PlantHeight', 'FieldWeight', 'BareCobWeight', 'DaysToFlowering')

## count up all non-deregressed accessions that were originally grown in (non-regressed) field trials
blups_std_acc_filter = blups_std %>% filter(Trial_Classification != 'Stress')

trial_trait_table = setNames(data.frame(matrix(ncol = 6, nrow = 23),
                                        row.names = blups_std_acc_filter$Experimento %>% unique()),
                             traits)
for(trait in traits) {
  for(experiment in blups_std_acc_filter$Experimento %>% unique()) {
    num_acc = blups_std_acc_filter %>% filter(Trait == trait, Experimento == experiment) %>% pull(SampleID) %>% unique() %>% length()
    trial_trait_table[experiment, trait] = num_acc
  }
}

write.csv(trial_trait_table, here(analyses_dir, 'Tables', 'accessions_per_trial_per_trait_original.csv'), quote = F)

## count up all deregressed accessions that were grown in (non-regressed) field trials AND had climate data that were used in our analyses
blups_deregressed_acc_filter = blups_deregressed %>% filter(Trial_Classification != 'Stress' & Trait != 'GrainWeightPerHectare') %>% filter(SampleID %in% finalMat$SampleID)

# trial_trait_table_deregressed = setNames(data.frame(matrix(ncol = 6, nrow = 23), 
#                                         row.names = blups_deregressed_acc_filter$Experimento %>% unique()),
#                              traits)
# for(trait in traits) {
#   for(experiment in blups_deregressed_acc_filter$Experimento %>% unique()) {
#     num_acc = blups_deregressed_acc_filter %>% filter(Trait == trait, Experimento == experiment) %>% pull(SampleID) %>% unique() %>% length()
#     trial_trait_table_deregressed[experiment, trait] = num_acc
#   }
# }
# write.csv(trial_trait_table_deregressed, here(analyses_dir, 'Tables', 'accessions_per_trial_per_trait_deregressed.csv'))

accessions_per_trial = blups_deregressed_acc_filter %>% group_by(Experimento) %>% summarize(n_unique = n_distinct(SampleID))

## how many accessions in each trial:tester combo?
accessions_per_tester_per_trial = table(blups_deregressed_acc_filter$Experimento, blups_deregressed_acc_filter$Tester)
write.csv(accessions_per_tester_per_trial, here(analyses_dir, 'Tables', 'accessions_per_trial_per_tester.csv'), quote = F)
accessions_per_trial_per_trait = table(blups_deregressed_acc_filter$Experimento, blups_deregressed_acc_filter$Trait)
write.csv(accessions_per_trial_per_trait, here(analyses_dir, 'Tables', 'accessions_per_trial_per_trait.csv'), quote = F)

## how many testers are represented in each trial:trait combo?
testers_in_trial = blups_deregressed_acc_filter %>% 
  group_by(Trait, Experimento) %>% 
  summarise(tester_count = n_distinct(Tester))
testers_in_trial = pivot_wider(testers_in_trial, names_from = Trait, values_from = tester_count, values_fill = list(tester_count = 0))
write.csv(testers_in_trial, here(analyses_dir, 'Tables', 'testers_in_trial.csv'), quote = F)

## how many trials are represented for each tester:trait combo?
trials_per_tester = blups_deregressed_acc_filter %>% 
  group_by(Trait, Tester) %>% 
  summarise(trial_count = n_distinct(Experimento))
trials_per_tester = pivot_wider(trials_per_tester, names_from = Trait, values_from = trial_count, values_fill = list(trial_count = 0))
write.csv(trials_per_tester, here(analyses_dir, 'Tables', 'trials_per_tester.csv'), quote = F)


# library(xlsx)
# write.xlsx(trial_trait_table_deregressed, file = 'SeeD_trial_information.xlsx', sheetName = 'accessions in trial per trait')
# write.xlsx(testers_in_trial, file = 'SeeD_trial_information.xlsx', sheetName = 'testers in trial')
# write.xlsx(trials_per_tester, file = 'SeeD_trial_information.xlsx', sheetName = 'trials per tester')
# write.xlsx(table(filtered_blups_deregressed$Experimento, filtered_blups_deregressed$Tester), 
#            file = 'SeeD_trial_information.xlsx', sheetName = 'accessions in trial per tester')

####################################################################################
## Just print out samples of matched SNPs and be done with this farce
####################################################################################
gea_matching_samples = as.data.frame(replicate(20, compare_GWAS_effects(gea_results, top_hits_clumped, gea = T) %>% filter(matching == T) %>% pull(SNP)))
write.table(gea_matching_samples, here(genetic_data_dir, 'gea_matching_sampled_SNPs.txt'), col.names = F, row.names = F)

## print out GEA SNPs under 1e-2 p-value for plink clumping and LD analysis
SNP_LD = gea_results %>% filter(P < 1e-2) %>% select(SNP)
write.csv(SNP_LD, 'Analyses/GEA_output/multivariate_results/GEA_SNP_1e2_forplink.txt', row.names = F, quote = F)

plink_analysis = fread('Analyses/GEA_output/multivariate_results/clumped/plink/LD-analysis.ld')

manual_manhattan(gea_results, plink_analysis$SNP_B)

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


write.csv(gea_beta_hats %>% filter(P < 1e-5), here(analyses_dir, 'Tables/envGWAS_significant_SNPs_effectsizes.csv'))

write_pheno_fx = function(pheno_process, pheno_name) {
  pheno_process = pheno_process %>% 
    filter(SNP %in% (gea_beta_hats %>% filter(P < 1e-5) %>% pull(SNP)))
  write.csv(pheno_process, here(analyses_dir, sprintf('Tables/significant_SNPs_%s_effectsizes.csv', pheno_name)))
}
purrr::map2(list(bcw_process, gwph_process, fw_process, ph_process, dtf_process, asi_process), 
                                    list('bare cob weight', 'grain weight per hectare', 'field weight', 'plant height', 'days to flowering', 'ASI'), 
                           write_pheno_fx)



####################################################################################
## phenotypic prediction
####################################################################################

## make tables for prediction ability values

## doesn't include stress trials or testers with low sample sizes
all_traits_prediction_ability_deregressed_nostress_filterTrial = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
  cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)), trait_i)
})))

# mex = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
#   cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_mex_nostress_filterTrial/modelPrediction_nostress_%s_results_resid.csv', trait_i)), trait_i)
# })))

# all_traits_prediction_ability_deregressed = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
#   cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups/modelPrediction_%s_results_resid.csv', trait_i)), trait_i)
# })))
# 
# all_traits_prediction_ability_deregressed_nostress = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
#   cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_nostress/modelPrediction_nostress_%s_results_resid.csv', trait_i)), trait_i)
# })))
# 
# all_traits_prediction_ability_deregressed_filterTrial = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
#   cbind(fread(sprintf('Analyses/PhenotypicPrediction/deregressed_blups_filterTrial/modelPrediction_%s_results_resid.csv', trait_i)), trait_i)
# })))
# 
# all_traits_prediction_ability_std = data.frame(data.table::rbindlist(lapply(traits[-c(1)], function(trait_i){
#   cbind(fread(sprintf('Analyses/PhenotypicPrediction/std/modelPrediction_%s_results_resid.csv', trait_i)), trait_i)
# })))
write.csv(all_traits_prediction_ability_deregressed_nostress_filterTrial, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_deregressed_nostress_filterTrial.csv'))
# write.csv(all_traits_prediction_ability_deregressed, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_deregressed.csv'))
# write.csv(all_traits_prediction_ability_std, here(analyses_dir, 'Tables', 'all_traits_prediction_ability_std.csv'))

## analyze mean prediction ability across trials
# all_traits_prediction_ability_deregressed_melt = pivot_longer(all_traits_prediction_ability_deregressed,
#                                                               cols = c("lm_PC5", 'lm_matching_SNPs',
#                                                                        'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
#                                                                        'lm_PC5_env', 'rf_env', 'rf_env_PC5'),
#                                                               names_to = 'model',
#                                                               values_to = 'r')


all_traits_prediction_ability_deregressed_melt = pivot_longer(all_traits_prediction_ability_deregressed_nostress_filterTrial, 
                                                              cols = c("lm_PC5", 'lm_matching_SNPs', 'lm_env',
                                                                       'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
                                                                       'lm_PC5_env', 'rf_env', 'rf_env_PC5'),
                                                              names_to = 'model',
                                                              values_to = 'r')
## just for mexican accessions
# all_traits_prediction_ability_deregressed_melt = pivot_longer(mex, 
#                                                               cols = c("lm_PC5", 'lm_matching_SNPs', 'lm_env',
#                                                                        'lm_enriched_SNPs', 'lm_all_SNPs', 'lm_all_SNPs_env',
#                                                                        'lm_PC5_env', 'rf_env', 'rf_env_PC5'),
#                                                               names_to = 'model',
#                                                               values_to = 'r')



## obtain prediction means for each trial
all_traits_prediction_ability_deregressed_trial = all_traits_prediction_ability_deregressed_nostress_filterTrial %>% 
  dplyr::select(-c(V1, n_train, n_test, n_total)) %>% 
  group_by(trait_i, Experimento) %>% 
  summarise_at(vars(-c(fold)),mean)
write.csv(all_traits_prediction_ability_deregressed_trial, here(analyses_dir, 'Tables', 'all_traits_mean_prediction_ability_trial.csv'), quote = F)

## take prediction means for entire trait across all trials
all_traits_prediction_ability_deregressed_trait = all_traits_prediction_ability_deregressed_nostress_filterTrial %>% 
  dplyr::select(-c(V1, n_train, n_test, n_total)) %>% 
  group_by(trait_i) %>% 
  summarise_at(vars(-c(fold, Experimento)),mean)
write.csv(all_traits_prediction_ability_deregressed_trait, here(analyses_dir, 'Tables', 'all_traits_mean_prediction_ability_trait.csv'), quote = F)

all_traits_prediction_ability_deregressed_trial %>% 
  mutate(climate_difference = lm_all_SNPs_env - lm_all_SNPs,
         climate_better = sapply(climate_difference, function(x) ifelse(abs(x) < 0.01, 0, sign(x)))) %>% 
  filter(trait_i == 'FieldWeight') %>% 
  dplyr::select(c(Experimento, climate_difference, climate_better))

wilcox.test(all_traits_prediction_ability_deregressed_trial %>% filter(trait_i == 'FieldWeight') %>% pull(lm_all_SNPs),
            all_traits_prediction_ability_deregressed_trial %>% filter(trait_i == 'FieldWeight') %>% pull(lm_all_SNPs_env),
            paired = T, alternative = 'two.sided')

## run significance tests between models
## for all traits
{
  all_traits_deregressed_melt = all_traits_prediction_ability_deregressed_melt %>% 
    filter(trait_i %in% c('PlantHeight')) %>% 
    filter(model %in% c('lm_PC5', 'lm_PC5_env',
                        'lm_matching_SNPs', 'lm_enriched_SNPs', 
                        'lm_all_SNPs', 'lm_all_SNPs_env'))
  prediction_comparison_model = lmer(r ~ Experimento + model + (1|model:Experimento) + (1|Experimento:fold) , 
                                     data = all_traits_deregressed_melt, weights = n_total)
  anova(prediction_comparison_model, ddf = 'K')
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1-0.05/1)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/1, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value)
  prediction_comparison_effects_summary
}
## just for field weight
{
  fieldweight_deregressed_melt = all_traits_prediction_ability_deregressed_melt %>% 
    filter(trait_i == 'FieldWeight') %>% 
    filter(model %in% c('lm_PC5', 'lm_PC5_env',
                        'lm_matching_SNPs', 'lm_enriched_SNPs', 
                        'lm_all_SNPs', 'lm_all_SNPs_env'))
  prediction_comparison_model = lmer(r ~ Experimento + (1|fold) + (1|model:fold) + (1|Experimento:fold) + model:Experimento + model, 
                                   data = fieldweight_deregressed_melt, weights = n_total)
  model_anova = anova(prediction_comparison_model, ddf = 'K')
  combined_pvalue = mean(model_anova$`Pr(>F)`[2:3])
  print(model_anova)
  print(combined_pvalue)
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1-0.05/1)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/1, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value)
  prediction_comparison_effects_summary
}

## looking at interaction effect - model:trial is fixed
{
  fieldweight_deregressed_melt = all_traits_prediction_ability_deregressed_melt %>% 
    filter(trait_i == 'FieldWeight') %>% 
    filter(model %in% c('lm_PC5', 'lm_PC5_env',
                        'lm_matching_SNPs', 'lm_enriched_SNPs', 
                        'lm_all_SNPs', 'lm_all_SNPs_env'))
  n_experiments = fieldweight_deregressed_melt %>% pull(Experimento) %>% unique() %>% length()
  prediction_comparison_model = lmer(r ~ Experimento + model + model:Experimento + (1|Experimento:fold) , 
                                     data = fieldweight_deregressed_melt, weights = n_total)
  model_anova = anova(prediction_comparison_model, ddf = 'K')
  print(model_anova)
  emmip(prediction_comparison_model, model ~ Experimento, CIs = F)
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', by = 'Experimento', infer = c(T, T), level = 1-0.05/6)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise', name = 'model_effect')
  regrouped_model_effects = update(prediction_comparison_effects,by = 'model_effect')
  experimento_model_comparison_effects = contrast(regrouped_model_effects, 'pairwise')
  experimento_model_comparison_effects_summary = summary(experimento_model_comparison_effects,infer = T,level = 1-0.05/n_experiments)
  # prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/n_experiments, infer = c(T, T))
  # prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value)
  # prediction_comparison_effects_summary
  experimento_model_comparison_effects_summary$p.value = pmin(1,experimento_model_comparison_effects_summary$p.value * n_experiments)
}

## looking at main effect - make trial:model random
{
  fieldweight_deregressed_melt = all_traits_prediction_ability_deregressed_melt %>% 
    filter(trait_i == 'FieldWeight') %>% 
    filter(model %in% c('lm_PC5', 'lm_PC5_env',
                        'lm_matching_SNPs', 'lm_enriched_SNPs', 
                        'lm_all_SNPs', 'lm_all_SNPs_env'))
  prediction_comparison_model = lmer(r ~ Experimento + model + (1|model:Experimento) + (1|Experimento:fold) , 
                                     data = fieldweight_deregressed_melt, weights = n_total)
  model_anova = anova(prediction_comparison_model, ddf = 'K')
  print(model_anova)
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1-0.05/1)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/1, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value)
  prediction_comparison_effects_summary
}
## for purely environmental
{
  fieldweight_deregressed_env_melt = all_traits_prediction_ability_deregressed_melt %>% 
    filter(trait_i == 'PlantHeight') %>% 
    filter(model %in% c('lm_PC5', 'lm_env',
                        'rf_env', 'rf_env_PC5'))
  # prediction_comparison_model = lmer(r ~ Experimento + (1|fold) + (1|model:fold) + (1|Experimento:fold) + model:Experimento + model, 
                                     # data = fieldweight_deregressed_env_melt)
  prediction_comparison_model = lmer(r ~ Experimento + model + (1|model:Experimento) + (1|Experimento:fold) , 
                                     data = fieldweight_deregressed_env_melt, weights = n_total)
  anova(prediction_comparison_model, ddf = 'K')
  prediction_comparison_means = emmeans(prediction_comparison_model, specs = 'model', infer = c(T, T), level = 1-0.05/1)
  cld(prediction_comparison_means)
  prediction_comparison_effects = contrast(prediction_comparison_means, 'pairwise')
  prediction_comparison_effects_summary = summary(prediction_comparison_effects, level = 1-0.05/1, infer = c(T, T))
  prediction_comparison_effects_summary$p.value = pmin(1,prediction_comparison_effects_summary$p.value)
  prediction_comparison_effects_summary
}

all_traits_prediction_ability_deregressed_melt %>% 
  group_by(trait_i, model) %>% 
  summarise(count = n(), mean_r = mean(r, na.rm = T), min_r = min(r), max_r= max(r)) %>% 
  filter(trait_i == 'FieldWeight') %>% 
  print(n = 1000)

t.test(all_traits_prediction_ability_deregressed_trial %>% filter(trait_i == 'BareCobWeight') %>% pull(rf_env), alternative = 'greater')
