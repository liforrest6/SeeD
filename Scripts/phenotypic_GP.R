####################################################################################
# How well does random forest work for predicting phenotypic values?
# Do environmental values work better than PCs?  Do testers determine all?
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################
library(data.table)
library(dplyr)
library(randomForest)
library(foreach)
library(doParallel)
library(stringr)
library(ggplot2)
library(tidyverse)
library(rrBLUP)
# library(lme4)
library(ModelMetrics)
registerDoParallel(4)

genetic_data_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data'
phenotype_data_dir = '/group/runciegrp2/Projects/SeeD/Phenotype_data'
env_data_dir = '/group/runciegrp2/Projects/SeeD/Env_data'
analyses_dir = '/group/runciegrp2/Projects/SeeD/Analyses'


# run = as.numeric(commandArgs(t=T)[1])
# run = as.numeric(strsplit(as.character(run),'')[[1]])
run = commandArgs(t=T)
# rep = run[1]
# dataset = run[2]
transform = run[1]
traitN = as.numeric(run[2])
ntree = as.numeric(run[3])
output = run[4]
nostress = run[5]
mex = run[6]


if(is.na(transform)) transform = 1
if(is.na(traitN)) traitN = 4
if(is.na(ntree)) ntree = 100
if(is.na(output)) output = 'deregressed_blups'
if(is.na(nostress)) nostress = ''
if(is.na(mex)) mex = ''


filterTrial = '_filterTrial'
print(output)
print(nostress)
print(filterTrial)

# transform = 1
# traitN = 4

K = fread(file.path(genetic_data_dir, 'Imputed_V4', 'K_allChr.csv'),data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
sample_to_geneticData = fread(file.path(genetic_data_dir, 'Imputed_V4', 'selected_genotypeIDs.csv'),data.table=F)

# data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
# trait = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")[traitN]
trait = c("DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","BareCobWeight","PlantHeight")[traitN]
# data$DNAID = sample_to_geneticData$V1[match(data$SampleID,sample_to_geneticData$Sample)]
# 
# data = subset(data,Trait == trait & !is.na(DNAID))

finalMat = read.csv(file.path(env_data_dir, 'GEA-climate-invnormtransformed.csv'))

# finalMat_race = read.csv(file.path(env_data_dir, 'GEA-race-classifications.csv'))
# finalMat = cbind(finalMat, PrimaryRace=finalMat_race$PrimaryRace)

#--------------------------------------------------------------------------------------------------------------------------------------------#
## test with only central mexican accessions
if(mex != '') {
  output = paste0(output, mex)
  filt_df = read.csv(file.path(env_data_dir, 'central_mexican_accessions.csv'))
  finalMat = finalMat %>% filter(Unique.ID %in% filt_df$Unique.ID.Mex)
}

just_clim = finalMat[-c(1)]
if(output == 'deregressed_blups') {
  blup_phenotypes = read.csv(file.path(phenotype_data_dir, 'blups_deregressed.csv'))
  } else {
    blup_phenotypes = read.csv(file.path(phenotype_data_dir, 'blups_std.csv'))
  }

## make output directories for plots and analyses
dir.create(file.path(sprintf('/group/runciegrp2/Projects/SeeD/Plots/phenotypic-consequence/%s', paste0(output, nostress, filterTrial))), showWarnings = F)
dir.create(file.path(analyses_dir, sprintf('PhenotypicPrediction/%s', paste0(output, nostress, filterTrial))), showWarnings = F)


finalMat$SampleID = str_split_i(finalMat$Unique.ID, ':', i = 1)
# envMat = separate(data = finalMat, col = Unique.ID, into = c('SampleID', 'Unique.ID'), sep = ':', remove = F)
#--------------------------------------------------------------------------------------------------------------------------------------------#

phenotypes_env = blup_phenotypes %>% filter(SampleID %in% finalMat$SampleID) %>% merge(finalMat, by = 'SampleID')
data = phenotypes_env

enriched_design = read.csv(file.path(genetic_data_dir, 'Imputed_V4', 'genotypes_by_chromosome', 'GEA_clumped_SNPs_genotype.csv'), row.names = 1)
matching_design = read.csv(file.path(genetic_data_dir, 'Imputed_V4', 'genotypes_by_chromosome', 'GEA_matching_SNPs_genotype.csv'), row.names = 1)
# rownames(enriched_design) = str_split_i(rownames(enriched_design), ':', i = 1)

print(traitN)
print(trait)

# remove trials with stress
if(nostress == '_nostress') {
  data = subset(data,Trial_Classification != 'Stress')
}

## this is cross-fold with OOB testing for each tester (test cross-tester predictive ability)
tester_results = foreach(exp = unique(data$Experimento),.combine = bind_rows, .errorhandling = 'remove') %dopar% {
  # exp = 'm2012BCH'
  # tester = 'CML269/CML264'

  data_env = subset(data,Experimento == exp & Trait == trait)

  ## remove accessions that are of CML244/CML349 if in thesetwo trials due to low sample size
  if(filterTrial == '_filterTrial') {
    if(exp == 'm2011BMO' | exp == 'm2011BTLNN') {
      data_env = subset(data_env, Tester != 'CML244/CML349')
    }
    if(exp == 'm2012BAL') {
      data_env = subset(data_env, Tester != 'CML451/CML486')
    }
  }

  residuals = resid(lm(Value ~ Tester, data = data_env))
  data_env$residuals = unname(residuals)
  print(exp)

  ## this does INT transform but everything should already have been done
  # if(transform == 2) {
  #   x = data_env$Value
  #   data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  # }

  ## calculate PCs from subset of accessions
  K_env = K[data_env$Unique.ID,data_env$Unique.ID]
  sK = svd(K_env)

  ## subset enriched GEA SNP genotypes by SampleID
  enriched_design_byTester = enriched_design[data_env$Unique.ID, ]
  matching_design_byTester = matching_design[data_env$Unique.ID, ]


  ## withholding one tester per fold
  fold_result = foreach(tester = unique(data_env$Tester),.combine=bind_rows, .errorhandling = 'remove') %do% {

    ### generate folds by withholding one tester
    fold_indices = data_env$Tester == tester
    data_train = data_env[!data_env$Tester == tester,]
    data_test = data_env[data_env$Tester == tester,]
    y_train = data_train$residuals
    y_test = data_test$residuals

    ## add PCs to env matrices
    sK_train = sK$u[!fold_indices, 1:5]
    sK_test = sK$u[fold_indices, 1:5]

    ## create model matrix for rrBLUP
    X = model.matrix(~Tester,data_env)
    y = data_env$residuals
    ## make sure that we don't have linearly dependent columns
    dropcols = caret::findLinearCombos(X[!fold_indices,])
    X_tester = X
    if(length(dropcols$remove) > 0) {
      X_tester = X[,-dropcols$remove]
    }
    ## withhold for testing
    y[fold_indices] = NA

    print(tester)


    ## only tester model
    # lm_res_X = mixed.solve(y,X = X_tester)

    ## model tester + first 5 PCs (pop structure)
    lm_PC5_res_X = mixed.solve(y, X = cbind(X_tester,sK$u[,1:5]))
    ## PCs + climate data
    lm_PC5_env_X = mixed.solve(y, X = cbind(X_tester, sK$u[,1:5], data_env %>% 
                                              select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))

    ## PCs + GEA SNPs
    lm_enriched_SNPs_X = mixed.solve(y, X = cbind(X_tester, sK$u[,1:5]), Z = enriched_design_byTester)
    lm_enriched_SNPs_only_X = mixed.solve(y, X = cbind(X_tester), Z = enriched_design_byTester)
    lm_matching_SNPs_X = mixed.solve(y, X = cbind(X_tester, sK$u[,1:5]), Z = matching_design_byTester)

    ## all GBS data and GBS + climate
    lm_all_SNPs_X = mixed.solve(y, X = cbind(X_tester, sK$u[,1:5]), K = K_env)
    lm_all_SNPs_env_X = mixed.solve(y, X = cbind(X_tester, sK$u[,1:5], 
      data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)), K = K_env)

    ## lmer with race as random effect
    # lm_race_X = lmer(residuals ~ (1|PrimaryRace), data_train)

    ## lmm with only climate data
    lm_env_X = mixed.solve(y, X = cbind(X_tester, data_env %>% 
                                              select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))

    ## random forest model with just environmental data
    env_Forest = randomForest(x = cbind(data_train %>% 
                                            select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)), 
                                y = y_train,
                                ntree=ntree,
                                importance=TRUE)
    env_res_X = predict(env_Forest, newdata = data_test %>% select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))

    ## random forest model with  environmental data and PCs
    env_PC5_Forest = randomForest(x = cbind(data_train %>% 
                                            select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK_train), 
                                y = y_train,
                                ntree=ntree,
                                importance=TRUE)
    env_PC5_res_X = predict(env_PC5_Forest, newdata = cbind(data_test %>% select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK_test))


    print(sprintf('Made models for %s', tester))

    # lm_tester_results = cor(cbind(data_env$residuals, X_tester %*% lm_res_X$beta)[fold_indices,])[2]
    lm_PC5_results = cor(cbind(data_env$residuals,cbind(X_tester,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[fold_indices,])[2]
    lm_PC5_env_results = cor(cbind(data_env$residuals,
                  as.matrix(cbind(X_tester,sK$u[,1:5], data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                  %*% lm_PC5_env_X$beta)[fold_indices,])[2]

    lm_enriched_SNPs_results = cor(cbind(data_env$residuals, c(cbind(X_tester,sK$u[,1:5]) %*% lm_enriched_SNPs_X$beta + 
      as.matrix(enriched_design_byTester) %*% matrix(lm_enriched_SNPs_X$u, ncol = 1)))[fold_indices,])[2]
    lm_enriched_SNPs_only_results = cor(cbind(data_env$residuals, c(X_tester %*% lm_enriched_SNPs_X$beta + 
      as.matrix(enriched_design_byTester) %*% matrix(lm_enriched_SNPs_only_X$u, ncol = 1)))[fold_indices,])[2]
    lm_matching_SNPs_results = cor(cbind(data_env$residuals, c(cbind(X_tester,sK$u[,1:5]) %*% lm_matching_SNPs_X$beta + 
      as.matrix(matching_design_byTester) %*% matrix(lm_matching_SNPs_X$u, ncol = 1)))[fold_indices,])[2]

    lm_all_SNPs_results = cor(cbind(data_env$residuals,cbind(X_tester,sK$u[,1:5]) %*% lm_all_SNPs_X$beta + matrix(lm_all_SNPs_X$u, ncol = 1))[fold_indices,])[2]
    lm_all_SNPs_env_results = cor(cbind(data_env$residuals,
      as.matrix(cbind(X_tester,sK$u[,1:5], data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
      %*% lm_all_SNPs_env_X$beta + matrix(lm_all_SNPs_env_X$u, ncol = 1))[fold_indices,])[2]

    lm_env_results = cor(cbind(data_env$residuals,
                  as.matrix(cbind(X_tester,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                  %*% lm_env_X$beta)[fold_indices,])[2]
    rf_env_results = cor(y_test, env_res_X)
    rf_env_PC5_results = cor(y_test, env_PC5_res_X)
    # lm_race_results = cor(y_test,predict(lm_race_X,newdata=data_test,allow.new.levels=TRUE))

    tester_result = data.frame(Experimento = exp, fold = tester, 
      n_train = length(y_train), n_test = length(y_test), n_total = length(y_train) + length(y_test),
               lm_PC5 = lm_PC5_results,
               lm_PC5_env = lm_PC5_env_results,
               lm_all_SNPs = lm_all_SNPs_results,
               lm_all_SNPs_env = lm_all_SNPs_env_results,
               lm_enriched_SNPs = lm_enriched_SNPs_results,
               lm_enriched_SNPs_only = lm_enriched_SNPs_only_results,
               lm_matching_SNPs = lm_matching_SNPs_results,
               lm_env = lm_env_results,
               rf_env = rf_env_results,
               rf_env_PC5 = rf_env_PC5_results
               # lm_race = lm_race_results,
               )
    print(sprintf('Tested models for %s', tester))

    tester_result
    
  }
  fold_result
}

results_file = sprintf(file.path(analyses_dir, 'PhenotypicPrediction', paste0(output, nostress, filterTrial), 'modelPrediction_sjh%s_%s_results_resid.csv'), nostress,trait)
print(results_file)
write.csv(tester_results,file = results_file)



# # make plots


# results_long = pivot_longer(results, cols = c(rf_tester, rf_PC5, rf_PC10, rf_env, rf_envPC, lm_K, lm_PC5, lm_env), names_to = 'model')
# plot_5fold = ggplot(results_long, aes(x = model, y = value, color = model)) +
#   geom_violin() +
#   facet_wrap(facets = 'Experimento') +
#   ylab('Predictive accuracy (r^2)') +
#   ggtitle('5-fold models accuracy improvement over only-tester model')
# png(sprintf('/group/runciegrp2/Projects/SeeD/Plots/RF/%s/RF_5fold_%s_resid%s.png', output, trait, mex), width = 655, height = 655)
# print(plot_5fold)
# dev.off()


# tester_results_long = pivot_longer(tester_results, cols = c(lm_PC5, lm_env, rf_env, lm_all_SNPs, lm_enriched_SNPs, lm_matching_SNPs), names_to = 'model')
# plot_byTester = ggplot(tester_results_long, aes(x = model, y = value, color = model)) +
#   geom_jitter() +
#   facet_wrap(facets = 'Experimento') +
#   ylab('Pearson correlation (r)') +
#   theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
#   ggtitle(sprintf('Predictive accuracy for %s model prediction within tester', trait))

# png(sprintf('/group/runciegrp2/Projects/SeeD/Plots/phenotypic-consequence/%s/modelPrediction_byTester_%s_resid%s.png', output, trait, mex), width = 655, height = 655)
# print(plot_byTester)
# dev.off()
