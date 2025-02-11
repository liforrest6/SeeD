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
mex = run[5]


if(is.na(transform)) transform = 1
if(is.na(traitN)) traitN = 4
if(is.na(ntree)) traitN = 1000
if(is.na(output)) output = 'resid'
if(is.na(mex)) mex = ''

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

## test with only central mexican accessions
if(mex != '') {
  output = paste0(output, mex)
  filt_df = read.csv(file.path(env_data_dir, 'central_mexican_accessions.csv'))
  finalMat = finalMat %>% filter(Unique.ID %in% filt_df$Unique.ID.Mex)
}

just_clim = finalMat[-c(1)]
blup_phenotypes = read.csv(file.path(phenotype_data_dir, 'blups_std.csv'))
finalMat$SampleID = str_split_i(finalMat$Unique.ID, ':', i = 1)
# envMat = separate(data = finalMat, col = Unique.ID, into = c('SampleID', 'Unique.ID'), sep = ':', remove = F)

phenotypes_env = blup_phenotypes %>% filter(SampleID %in% finalMat$SampleID) %>% merge(finalMat, by = 'SampleID')
data = phenotypes_env

print(traitN)
print(trait)

## this is 5-fold cross-validation within each trial
results = foreach(exp = unique(data$Experimento),.combine = bind_rows, .errorhandling = 'remove') %dopar% {
  
  # for(exp in unique(data$Experimento)) {
  data_env = subset(data,Experimento == exp & Trait == trait)
  residuals = resid(lm(Value ~ Tester, data = data_env))
  data_env$residuals = unname(residuals)
  # if(dim(data_env)[1] == 0) {
  #   next
  # }
  ## this does INT transform but everything should already have been done
  # if(transform == 2) {
  #   x = data_env$Value
  #   data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  # }

  K_env = K[data_env$Unique.ID,data_env$Unique.ID]
  sK = svd(K_env)

  # X = model.matrix(~Tester,data_env)
  # X_PC = cbind(X,sK$u[,1:10])

  folds = caret::createFolds(data_env$Value,k=5)

  env_results = foreach(i_fold = seq_along(folds),.combine=bind_rows) %do%   {
    
    fold_indices = rownames(data_env)[folds[[i_fold]]]
    data_train = data_env[!rownames(data_env) %in% fold_indices,]
    data_test = data_env[rownames(data_env) %in% fold_indices,]
    y_train = data_train$residuals
    y_test = data_test$residuals
    
    testerForest = randomForest(x = data_train %>% select(Tester), y = y_train,
                                ntree=ntree,
                                importance=TRUE)
    tester_res_X = predict(testerForest, newdata = data_test %>% select(Tester))
    
    # residuals = resid(lm(Value ~ Tester, data = data_env))
    # resid_train = data.frame(residuals[!names(residuals) %in% fold_indices])
    # resid_test = data.frame(residuals[names(residuals) %in% fold_indices])
    # colnames(resid_train) = 'residual'
    # colnames(resid_test) = 'residual'

    # testerResidForest = randomForest(x = resid_train, y = y_train,
    #                             ntree=ntree,
    #                             importance=TRUE)
    # testerResid_res_X = predict(testerResidForest, newdata = resid_test)
    
    PC5_Forest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:5]), y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    PC5_res_X = predict(PC5_Forest, newdata = cbind(data_test %>% select(Tester), sK$u[folds[[i_fold]],1:5]))
    
    PC10_Forest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:10]), y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    PC10_res_X = predict(PC10_Forest, newdata = cbind(data_test %>% select(Tester), sK$u[folds[[i_fold]],1:10]))
    
    # PC5_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:5]), y = y_train,
    #                           ntree=ntree,
    #                           importance=TRUE)
    # PC5_tester_res_X = predict(PC5_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[folds[[i_fold]],1:5]))
    
    # PC10_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:10]), y = y_train,
    #                            ntree=ntree,
    #                            importance=TRUE)
    # PC10_tester_res_X = predict(PC10_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[folds[[i_fold]],1:10]))
    
    env_Forest = randomForest(x = data_train %>% 
                                          select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), 
                              y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    env_res_X = predict(env_Forest, newdata = data_test %>% 
                                                      select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))
    
    envPC_Forest = randomForest(x = cbind(data_train %>% 
                                          select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[-folds[[i_fold]],1:5]), 
                              y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    envPC_res_X = predict(envPC_Forest, newdata = cbind(data_test %>% 
                                                      select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[folds[[i_fold]],1:5]))
    
    # allForest = randomForest(x = cbind(data_train %>% 
    #                                             select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[-folds[[i_fold]],1:5]), 
    #                                 y = y_train,
    #                                 ntree=ntree,
    #                                 importance=TRUE)
    # all_res_X = predict(allForest, newdata = cbind(data_test %>% 
    #                                            select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[folds[[i_fold]],1:5]))
    
    X = model.matrix(~Tester,data_env)
    y = data_env$residuals
    y[folds[[i_fold]]] = NA

    ## drop columns that are not full rank (can happen when only one accession in a tester, etc)
    dropcols = caret::findLinearCombos(X[!rownames(data_env) %in% fold_indices,])
    if(length(dropcols$remove) > 0) {
      X = X[,-dropcols$remove]
    }

    lm_K_res_X = mixed.solve(y,X = X,K=K_env)
    lm_env_res_X = mixed.solve(y, X = cbind(X, data_env %>% 
                                              select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))
    lm_PC5_res_X = mixed.solve(y, X = cbind(X,sK$u[,1:5]))

    print(exp)
    baseline = cor(y_test, tester_res_X)
    data.frame(Experimento = exp, fold = i_fold,
               rf_tester = cor(y_test, tester_res_X) - baseline,
               # testerResid = cor(y_test, testerResid_res_X) - baseline,
               rf_PC5 = cor(y_test, PC5_res_X) - baseline,
               rf_PC10 = cor(y_test, PC10_res_X) - baseline,
               # tester_PC5 = cor(y_test, PC5_tester_res_X)^2 - baseline,
               # tester_PC10 = cor(y_test, PC10_tester_res_X)^2 - baseline,
               rf_env = cor(y_test, env_res_X) - baseline,
               rf_envPC = cor(y_test, envPC_res_X) - baseline,
               # allVar = cor(y_test, all_res_X)^2 - baseline
               lm_K = cor(cbind(data_env$residuals,c(X %*% lm_K_res_X$beta) + lm_K_res_X$u)[folds[[i_fold]],])[2] - baseline,
               lm_PC5 = cor(cbind(data_env$residuals,cbind(X,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[folds[[i_fold]],])[2] - baseline,
               lm_env = cor(cbind(data_env$residuals,
                                  as.matrix(cbind(X,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                                  %*% lm_env_res_X$beta)[folds[[i_fold]],])[2] - baseline
               )
    
  }
  env_results
}


write.csv(results,file = sprintf(file.path(analyses_dir, 'RandomForest', output, 'RF_5fold_%s_results_resid%s.csv'), trait, mex))

print('start tester')
## this is cross-fold with OOB testing for each tester (test cross-tester predictive ability)
tester_results = foreach(exp = unique(data$Experimento),.combine = bind_rows, .errorhandling = 'remove') %dopar% {
  # exp = 'm2012BCH'
  # tester = 'CML269/CML264'
  data_env = subset(data,Experimento == exp & Trait == trait)
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


  ## withholding one tester per fold
  fold_result = foreach(tester = unique(data_env$Tester),.combine=bind_rows) %do% {

    fold_indices = data_env$Tester == tester
    data_train = data_env[!data_env$Tester == tester,]
    data_test = data_env[data_env$Tester == tester,]
    y_train = data_train$residuals
    y_test = data_test$residuals
    
    testerForest = randomForest(x = data_train %>% select(Tester), y = y_train,
                                ntree=ntree,
                                importance=TRUE)
    tester_res_X = predict(testerForest, newdata = data_test %>% select(Tester))
    
    PC5_Forest = randomForest(x = sK$u[!fold_indices,1:5], y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    PC5_res_X = predict(PC5_Forest, newdata = sK$u[fold_indices,1:5])
    
    PC10_Forest = randomForest(x = sK$u[!fold_indices,1:10], y = y_train,
                               ntree=ntree,
                               importance=TRUE)
    PC10_res_X = predict(PC10_Forest, newdata = sK$u[fold_indices,1:10])
    
    # PC5_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[!fold_indices,1:5]), y = y_train,
    #                                 ntree=ntree,
    #                                 importance=TRUE)
    # PC5_tester_res_X = predict(PC5_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[fold_indices,1:5]))

    # PC10_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[!fold_indices,1:10]), y = y_train,
    #                                  ntree=ntree,
    #                                  importance=TRUE)
    # PC10_tester_res_X = predict(PC10_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[fold_indices,1:10]))
    
    env_Forest = randomForest(x = data_train %>% 
                                select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), 
                              y = y_train,
                              ntree=ntree,
                              importance=TRUE)
    env_res_X = predict(env_Forest, newdata = data_test %>% 
                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))
    
    envPC_Forest = randomForest(x = cbind(data_train %>% 
                                            select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[!fold_indices,1:5]), 
                                y = y_train,
                                ntree=ntree,
                                importance=TRUE)
    envPC_res_X = predict(envPC_Forest, newdata = cbind(data_test %>% 
                                                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[fold_indices,1:5]))
    
    # allForest = randomForest(x = cbind(data_train %>% 
    #                                      select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[!fold_indices,1:5]), 
    #                          y = y_train,
    #                          ntree=ntree,
    #                          importance=TRUE)
    # all_res_X = predict(allForest, newdata = cbind(data_test %>% 
    #                                                  select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[fold_indices,1:5]))
    
    X = model.matrix(~Tester,data_env)
    y = data_env$residuals
    dropcols = caret::findLinearCombos(X[!fold_indices,])
    X_tester = X
    if(length(dropcols$remove) > 0) {
      X_tester = X[,-dropcols$remove]
    }
    
    
    y[fold_indices] = NA
    lm_K_res_X = mixed.solve(y,X = X_tester,K=K_env)
    lm_env_res_X = mixed.solve(y, X = cbind(X_tester, data_env %>% 
                                              select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))
    lm_PC5_res_X = mixed.solve(y, X = cbind(X_tester,sK$u[,1:5]))

    print(exp)
    # baseline = cor(y_test, tester_res_X)
    data.frame(Experimento = exp, fold = tester,
               # rf_tester = cor(y_test, tester_res_X),
               rf_PC5 = cor(y_test, PC5_res_X),
               rf_PC10 = cor(y_test, PC10_res_X),
               # tester_PC5 = cor(y_test, PC5_tester_res_X) - baseline,
               # tester_PC10 = cor(y_test, PC10_tester_res_X) - baseline,
               rf_env = cor(y_test, env_res_X),
               rf_envPC = cor(y_test, envPC_res_X),
               # allVar = cor(y_test, all_res_X) - baseline
               lm_K = cor(cbind(data_env$residuals,c(X_tester %*% lm_K_res_X$beta) + lm_K_res_X$u)[fold_indices,])[2],
               lm_PC5 = cor(cbind(data_env$residuals,cbind(X_tester,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[fold_indices,])[2],
               lm_env = cor(cbind(data_env$residuals,
                                  as.matrix(cbind(X_tester,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                                  %*% lm_env_res_X$beta)[fold_indices,])[2]
               )
    
  }
  fold_result
}

write.csv(tester_results,file = sprintf(file.path(analyses_dir, 'RandomForest', output, 'RF_byTester_%s_results_resid%s.csv'), trait, mex))



# # make plots


results_long = pivot_longer(results, cols = c(rf_tester, rf_PC5, rf_PC10, rf_env, rf_envPC, lm_K, lm_PC5, lm_env), names_to = 'model')
plot_5fold = ggplot(results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'Experimento') +
  ylab('Predictive accuracy (r^2)') +
  ggtitle('5-fold models accuracy improvement over only-tester model')
png(sprintf('/group/runciegrp2/Projects/SeeD/Plots/RF/%s/RF_5fold_%s_resid%s.png', output, trait, mex), width = 655, height = 655)
print(plot_5fold)
dev.off()


tester_results_long = pivot_longer(tester_results, cols = c(rf_PC5, rf_PC10, rf_env, rf_envPC, lm_K, lm_PC5, lm_env), names_to = 'model')
plot_byTester = ggplot(tester_results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'Experimento') +
  ylab('Predictive accuracy (r^2)') +
  ggtitle('Tester fold models raw predictive accuracy')

png(sprintf('/group/runciegrp2/Projects/SeeD/Plots/RF/%s/RF_byTester_%s_resid%s.png', output, trait, mex), width = 655, height = 655)
print(plot_byTester)
dev.off()
