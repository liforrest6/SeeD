####################################################################################
# Process results from multivariate GEA analyses
#
# Author: Forrest Li
# EnvGWAS analysis
####################################################################################

source(here::here('config.R'))
library(randomForest)

gea_results = read_GEA_results(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_results')
gea_beta_hats = read_GWAS_effects(here(analyses_dir, 'GEA_output', 'multivariate_results'), 'envGWAS_beta_hats')
gea_beta_hats = cbind(gea_results, gea_beta_hats)


top_hits_clumped = load_clumped()$SNP

qqman::qq(gea_results$P, main = 'qqplot, maf > 0.01')

qqman::qq(gea_results %>% filter(maf > 0.05) %>% pull(P), main = 'qqplot, maf > 0.05')

## read phenotypic joint GWAS results and effect sizes
{
  dtf_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'JointGWAS_output'), 'GWAS_results')
  dtf_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
  asi_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'ASI', 'JointGWAS_output'), 'GWAS_results')
  asi_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'ASI', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
  gwph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'JointGWAS_output'), 'GWAS_results')
  gwph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
  ph_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'JointGWAS_output'), 'GWAS_results')
  ph_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
  fw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'JointGWAS_output'), 'GWAS_results')
  fw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
  bcw_results = read_GWAS_results(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'JointGWAS_output'), 'GWAS_results')
  bcw_betas = read_GWAS_effects(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'JointGWAS_output'), 'GWAS_beta_hats_chr')
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
  dtf_SE = process_GWAS(dtf_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'DaysToFlowering', 'JointGWAS_output'), 'GWAS_SEs'))
  asi_SE = process_GWAS(asi_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'ASI', 'JointGWAS_output'), 'GWAS_SEs'))
  gwph_SE = process_GWAS(gwph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'GrainWeightPerHectareCorrected', 'JointGWAS_output'), 'GWAS_SEs'))
  ph_SE = process_GWAS(ph_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'PlantHeight', 'JointGWAS_output'), 'GWAS_SEs'))
  fw_SE = process_GWAS(fw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'FieldWeight', 'JointGWAS_output'), 'GWAS_SEs'))
  bcw_SE = process_GWAS(bcw_results, read_GWAS_SEs(here(analyses_dir, 'JointGWAS', 'BareCobWeight', 'JointGWAS_output'), 'GWAS_SEs'))
  
}


####################################################################################
# Regression testing
####################################################################################

gea_phenotype_fx = merge(gwph_process, gea_beta_hats, by = 'SNP')

summary(lm(Experimentom2011BAFNN.X ~ tmin..X + tmax..X + trange..X + precipTot..X + aridityMean..X + rhMean..X + elevation..X,
   data = gea_phenotype_fx))
summary(lm(Experimentom2011BCL.X ~ tmin..X + tmax..X + trange..X + precipTot..X + aridityMean..X + rhMean..X + elevation..X,
           data = gea_phenotype_fx))
summary(lm(Experimentom2011BMO.X ~ tmin..X + tmax..X + trange..X + precipTot..X + aridityMean..X + rhMean..X + elevation..X,
           data = gea_phenotype_fx))
summary(lm(Experimentom2012AOBRN.X ~ tmin..X + tmax..X + trange..X + precipTot..X + aridityMean..X + rhMean..X + elevation..X,
           data = gea_phenotype_fx))
summary(lm(Experimentom2012BEB.X ~ tmin..X + tmax..X + trange..X + precipTot..X + aridityMean..X + rhMean..X + elevation..X,
           data = gea_phenotype_fx))


####################################################################################
# Run random forest to test relative contribution of each environmental variable or tester 
####################################################################################
finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
just_clim = finalMat[-c(1)]
blup_phenotypes = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
envMat = separate(data = finalMat, col = Unique.ID, into = c('SampleID', 'Unique.ID'), sep = ':') %>% subset(select = -c(Unique.ID))

phenotypes_env = blup_phenotypes %>% filter(SampleID %in% envMat$SampleID) %>% merge(envMat, by = 'SampleID')


## GWPH
for(experimento in phenotypes_env %>% filter(Trait == 'GrainWeightPerHectareCorrected') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
               data = phenotypes_env %>% filter(Trait == 'GrainWeightPerHectareCorrected', Experimento == experimento), ntree=1000,
               keep.forest=FALSE, importance=TRUE)
  print(forest)
}

## PH
for(experimento in phenotypes_env %>% filter(Trait == 'PlantHeight') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
                        data = phenotypes_env %>% filter(Trait == 'PlantHeight', Experimento == experimento), ntree=1000,
                        keep.forest=FALSE, importance=TRUE)
  print(forest)
}

## ASI
for(experimento in phenotypes_env %>% filter(Trait == 'ASI') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
                        data = phenotypes_env %>% filter(Trait == 'ASI', Experimento == experimento), ntree=1000,
                        keep.forest=FALSE, importance=TRUE)
  print(forest)
}

## BCW
for(experimento in phenotypes_env %>% filter(Trait == 'BareCobWeight') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
                        data = phenotypes_env %>% filter(Trait == 'BareCobWeight', Experimento == experimento), ntree=1000,
                        keep.forest=FALSE, importance=TRUE)
  print(forest)
}

## DaysToFlowering
for(experimento in phenotypes_env %>% filter(Trait == 'DaysToFlowering') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
                        data = phenotypes_env %>% filter(Trait == 'DaysToFlowering', Experimento == experimento), ntree=1000,
                        keep.forest=FALSE, importance=TRUE)
  print(forest)
}

## FW
for(experimento in phenotypes_env %>% filter(Trait == 'FieldWeight') %>% pull(Experimento) %>% unique()){
  print(experimento)
  forest = randomForest(Value ~ tmin + tmax + trange + precipTot + aridityMean + rhMean + elevation + Tester, 
                        data = phenotypes_env %>% filter(Trait == 'FieldWeight', Experimento == experimento), ntree=1000,
                        keep.forest=FALSE, importance=TRUE)
  print(forest)
}


####################################################################################
## write a table for seeing how many accessions are in each trial for each trait
####################################################################################
phenotypes = c('GrainWeightPerHectareCorrected', 'PlantHeight', 'ASI', 'FieldWeight', 'BareCobWeight', 'DaysToFlowering')
trial_trait_table = setNames(data.frame(matrix(ncol = 6, nrow = 31), 
                                        row.names = blup_phenotypes$Experimento %>% unique()),
                             phenotypes)
for(trait in phenotypes) {
  for(experiment in blup_phenotypes$Experimento %>% unique()) {
    num_acc = phenotypes_env %>% filter(Trait == trait, Experimento == experiment) %>% pull(SampleID) %>% unique() %>% length()
    trial_trait_table[experiment, trait] = num_acc
  }
}
 
write.csv(trial_trait_table, here(phenotype_data_dir, 'Accessions_per_trial_per_trait.csv'))
          