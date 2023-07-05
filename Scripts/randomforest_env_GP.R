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
library(rrBLUP)
library(ModelMetrics)
registerDoParallel(4)

source(here::here('config.R'))


run = as.numeric(commandArgs(t=T)[1])
run = as.numeric(strsplit(as.character(run),'')[[1]])
# rep = run[1]
# dataset = run[2]
transform = run[1]
traitN = run[2]
ntree = run[3]

if(is.na(transform)) transform = 1
if(is.na(traitN)) traitN = 4

transform = 1
traitN = 4
ntree = 200

K = fread(here(genetic_data_dir, 'Imputed_V4', 'K_allChr.csv'),data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
sample_to_geneticData = fread(here(genetic_data_dir, 'Imputed_V4', 'selected_genotypeIDs.csv'),data.table=F)

# data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
trait = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")[traitN]
# data$DNAID = sample_to_geneticData$V1[match(data$SampleID,sample_to_geneticData$Sample)]
# 
# data = subset(data,Trait == trait & !is.na(DNAID))

finalMat = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))

filt_df = read.csv(here(env_data_dir, 'central_mexican_accessions.csv'))
# finalMat = finalMat %>% filter(Unique.ID %in% filt_df$Unique.ID.Mex)
just_clim = finalMat[-c(1)]
blup_phenotypes = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
finalMat$SampleID = str_split_i(finalMat$Unique.ID, ':', i = 1)

phenotypes_env = blup_phenotypes %>% filter(SampleID %in% finalMat$SampleID) %>% merge(finalMat, by = 'SampleID')
data = phenotypes_env

## this is 5-fold cross-validation within each trial
results = foreach(exp = unique(data$Experimento),.combine = bind_rows, .errorhandling = 'remove') %dopar% {
  exp = 'm2012BEB'

  data_env = subset(data,Experimento == exp & Trait == trait)
  residuals = resid(lm(Value ~ Tester, data = data_env))
  data_env$residuals = unname(residuals)


  K_env = K[data_env$Unique.ID,data_env$Unique.ID]
  sK = svd(K_env)

  folds = caret::createFolds(data_env$Value,k=5)

  env_results = foreach(i_fold = seq_along(folds),.combine=bind_rows) %do%   {
    
    fold_indices = rownames(data_env)[folds[[i_fold]]]
    data_train = data_env[!rownames(data_env) %in% fold_indices,]
    data_test = data_env[rownames(data_env) %in% fold_indices,]
    y_train = data_train$residuals
    y_test = data_test$residuals
    
    testerForest = randomForest(x = data_train %>% select(Tester), y = y_train,
                                ntree=200,
                                importance=TRUE)
    tester_res_X = predict(testerForest, newdata = data_test %>% select(Tester))
    # 
    # 
    # resid_train = data.frame(residuals[!names(residuals) %in% fold_indices])
    # resid_test = data.frame(residuals[names(residuals) %in% fold_indices])
    # colnames(resid_train) = 'residual'
    # colnames(resid_test) = 'residual'

    # testerResidForest = randomForest(x = resid_train, y = y_train,
    #                             ntree=200,
    #                             importance=TRUE)
    # testerResid_res_X = predict(testerResidForest, newdata = resid_test)
    
    PC5_Forest = randomForest(x = sK$u[-folds[[i_fold]],1:5], y = y_train,
                              ntree=200,
                              importance=TRUE)
    PC5_res_X = predict(PC5_Forest, newdata = sK$u[folds[[i_fold]],1:5])
    
    PC10_Forest = randomForest(x = sK$u[-folds[[i_fold]],1:10], y = y_train,
                              ntree=200,
                              importance=TRUE)
    PC10_res_X = predict(PC10_Forest, newdata = sK$u[folds[[i_fold]],1:10])
    
    # PC5_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:5]), y = y_train,
    #                           ntree=200,
    #                           importance=TRUE)
    # PC5_tester_res_X = predict(PC5_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[folds[[i_fold]],1:5]))
    # 
    # PC10_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[-folds[[i_fold]],1:10]), y = y_train,
    #                            ntree=200,
    #                            importance=TRUE)
    # PC10_tester_res_X = predict(PC10_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[folds[[i_fold]],1:10]))
    
    env_Forest = randomForest(x = data_train %>% 
                                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), 
                              y = y_train,
                              ntree=200,
                              importance=TRUE)
    env_res_X = predict(env_Forest, newdata = data_test %>% 
                                                      select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))
    
    envPC_Forest = randomForest(x = cbind(data_train %>% 
                                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[-folds[[i_fold]],1:5]), 
                              y = y_train,
                              ntree=200,
                              importance=TRUE)
    envPC_res_X = predict(envPC_Forest, newdata = cbind(data_test %>% 
                                                      select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[folds[[i_fold]],1:5]))
    
    # allForest = randomForest(x = cbind(data_train %>% 
    #                                             select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[-folds[[i_fold]],1:5]), 
    #                                 y = y_train,
    #                                 ntree=200,
    #                                 importance=TRUE)
    # all_res_X = predict(allForest, newdata = cbind(data_test %>% 
    #                                            select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[folds[[i_fold]],1:5]))
    # 
    
    X = model.matrix(~Tester,data_env)
    y = data_env$residuals
    y[folds[[i_fold]]] = NA

    dropcols = caret::findLinearCombos(X[!rownames(data_env) %in% fold_indices,])
    if(length(dropcols$remove) > 0) {
      X = X[,-dropcols$remove]
    }
    
    lm_K_res_X = mixed.solve(y,X = X,K=K_env)
    lm_env_res_X = mixed.solve(y, X = cbind(X, data_env %>% 
                                              select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))
    lm_PC5_res_X = mixed.solve(y, X = cbind(X,sK$u[,1:5]))

    
    baseline = cor(y_test, tester_res_X)^2
    data.frame(Experimento = exp, fold = i_fold,
               tester = cor(y_test, tester_res_X)^2 - baseline,
               # testerResid = cor(y_test, testerResid_res_X)^2 - baseline,
               PC5 = cor(y_test, PC5_res_X)^2 - baseline,
               PC10 = cor(y_test, PC10_res_X)^2 - baseline,
               # tester_PC5 = cor(y_test, PC5_tester_res_X)^2 - baseline,
               # tester_PC10 = cor(y_test, PC10_tester_res_X)^2 - baseline,
               env = cor(y_test, env_res_X)^2 - baseline,
               envPC = cor(y_test, envPC_res_X)^2 - baseline,
               lm_K = cor(cbind(data_env$residuals,c(X %*% lm_K_res_X$beta) + lm_K_res_X$u)[folds[[i_fold]],])[2] - baseline,
               lm_PC5 = cor(cbind(data_env$residuals,cbind(X,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[folds[[i_fold]],])[2] - baseline,
               lm_env = cor(cbind(data_env$residuals,
                                  as.matrix(cbind(X,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                                  %*% lm_env_res_X$beta)[folds[[i_fold]],])[2] - baseline
               # allVar = cor(y_test, all_res_X)^2 - baseline
               )
    
  }
  env_results
}


# write.csv(results,file = sprintf('single_trait_env_GP/GP_results_%s_%s.csv',trait,ifelse(transform == 1,'orig','IND')))

## this is cross-fold with OOB testing for each tester (test cross-tester predictive ability)
tester_results = foreach(exp = unique(data$Experimento),.combine = bind_rows, .errorhandling = 'remove') %dopar% {
  exp = 'm2012BCH'
  tester = 'CML269/CML264'
  data_env = subset(data,Experimento == exp & Trait == trait)
  residuals = resid(lm(Value ~ Tester, data = data_env))
  data_env$residuals = residuals
  ## this does INT transform but everything should already have been done
  # if(transform == 2) {
  #   x = data_env$Value
  #   data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  # }

  ## calculate PCs from subset of accessions
  K_env = K[data_env$Unique.ID,data_env$Unique.ID]
  sK = svd(K_env)
  
  # X = model.matrix(~Tester,data_env)
  # y = data_env$residuals
  # y[fold_indices] = NA

  ## withholding one tester per fold
  fold_result = foreach(tester = unique(data_env$Tester),.combine=bind_rows) %do% {
    # X = model.matrix(~Tester,data_env)
    
    fold_indices = data_env$Tester == tester
    data_train = data_env[!data_env$Tester == tester,]
    data_test = data_env[data_env$Tester == tester,]
    y_train = data_train$residuals
    y_test = data_test$residuals
    
    testerForest = randomForest(x = data_train %>% select(Tester), y = y_train,
                                ntree=200,
                                importance=TRUE)
    tester_res_X = predict(testerForest, newdata = data_test %>% select(Tester))
    
    PC5_Forest = randomForest(x = sK$u[!fold_indices,1:5], y = y_train,
                              ntree=200,
                              importance=TRUE)
    PC5_res_X = predict(PC5_Forest, newdata = sK$u[fold_indices,1:5])
    
    PC10_Forest = randomForest(x = sK$u[!fold_indices,1:10], y = y_train,
                               ntree=200,
                               importance=TRUE)
    PC10_res_X = predict(PC10_Forest, newdata = sK$u[fold_indices,1:10])
    
    # PC5_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[!fold_indices,1:5]), y = y_train,
    #                                 ntree=200,
    #                                 importance=TRUE)
    # PC5_tester_res_X = predict(PC5_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[fold_indices,1:5]))
    # 
    # PC10_testerForest = randomForest(x = cbind(data_train %>% select(Tester), sK$u[!fold_indices,1:10]), y = y_train,
    #                                  ntree=200,
    #                                  importance=TRUE)
    # PC10_tester_res_X = predict(PC10_testerForest, newdata = cbind(data_test %>% select(Tester),sK$u[fold_indices,1:10]))
    
    env_Forest = randomForest(x = data_train %>% 
                                select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), 
                              y = y_train,
                              ntree=200,
                              importance=TRUE)
    env_res_X = predict(env_Forest, newdata = data_test %>% 
                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))
    
    envPC_Forest = randomForest(x = cbind(data_train %>% 
                                            select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[!fold_indices,1:5]), 
                                y = y_train,
                                ntree=200,
                                importance=TRUE)
    envPC_res_X = predict(envPC_Forest, newdata = cbind(data_test %>% 
                                                          select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation),sK$u[fold_indices,1:5]))
    
    # allForest = randomForest(x = cbind(data_train %>% 
    #                                      select(Tester, tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation), sK$u[!fold_indices,1:5]), 
    #                          y = y_train,
    #                          ntree=200,
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
    
    # lm_K_res_X = mixed.solve(y,X = X,K=K_env)
    # lm_env_res_X = mixed.solve(y, X = cbind(X, data_env %>% 
    #                                           select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)))
    # lm_PC5_res_X = mixed.solve(y, X = cbind(X,sK$u[,1:5]))
    
    baseline = rmse(y_test, tester_res_X)
    # data.frame(Experimento = exp, fold = tester,
    #            rf_tester = cor(y_test, tester_res_X)^2,
    #            rf_PC5 = cor(y_test, PC5_res_X)^2 - baseline,
    #            rf_PC10 = cor(y_test, PC10_res_X)^2 - baseline,
    #            # tester_PC5 = cor(y_test, PC5_tester_res_X)^2 - baseline,
    #            # tester_PC10 = cor(y_test, PC10_tester_res_X)^2 - baseline,
    #            rf_env = cor(y_test, env_res_X)^2 - baseline,
    #            rf_envPC = cor(y_test, envPC_res_X)^2 - baseline,
    #            # allVar = cor(y_test, all_res_X)^2 - baseline
    #            lm_K = cor(cbind(data_env$residuals,c(X_tester %*% lm_K_res_X$beta) + lm_K_res_X$u)[fold_indices,])[2]^2 - baseline,
    #            lm_PC5 = cor(cbind(data_env$residuals,cbind(X_tester,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[fold_indices,])[2]^2 - baseline,
    #            lm_env = cor(cbind(data_env$residuals,
    #                               as.matrix(cbind(X_tester,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
    #                               %*% lm_env_res_X$beta)[fold_indices,])[2]^2 - baseline
    # )
    
    data.frame(Experimento = exp, fold = tester,
               rf_tester = rmse(y_test, tester_res_X),
               rf_PC5 = rmse(y_test, PC5_res_X),
               rf_PC10 = rmse(y_test, PC10_res_X),
               # tester_PC5 = cor(y_test, PC5_tester_res_X)^2 - baseline,
               # tester_PC10 = cor(y_test, PC10_tester_res_X)^2 - baseline,
               rf_env = rmse(y_test, env_res_X),
               rf_envPC = rmse(y_test, envPC_res_X),
               # allVar = cor(y_test, all_res_X)^2 - baseline
               lm_K = rmse(y_test, (c(X_tester %*% lm_K_res_X$beta) + lm_K_res_X$u)[fold_indices]),
               lm_PC5 = rmse(y_test, (cbind(X_tester,sK$u[,1:5]) %*% lm_PC5_res_X$beta)[fold_indices]) ,
               lm_env = rmse(y_test,
                                  (as.matrix(cbind(X_tester,data_env %>% select(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation))) 
                                  %*% lm_env_res_X$beta)[fold_indices])
    )
    
  }
  fold_result
}

# write.csv(results,file = sprintf('single_trait_env_GP/GP_tester_results_%s_%s.csv',trait,ifelse(transform == 1,'orig','INT')))



# # make plots
#
# # get H2 info
# traits = c('asi','dff','pcampo','pgrhah','polote','altpl','pgrha')
# traitNames = c('ASI','DaysToFlowering','FieldWeight','GrainWeightPerHectareCorrected','BareCobWeight','PlantHeight','GrainWeightPerHectare')
# H2s = as.data.frame(readxl::read_excel('../../Phenotype_data/Burgueno/MASAGRO-Maize-Data.xlsx',sheet=3))
# H2s$h2 = 1-H2s$`Average blup variance`/H2s$`Genetic variance`
# data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
# H2s$Experimento = unique(data$Experimento)[match(H2s$Experiment,toupper(unique(data$Experimento)))]
# H2s$TraitName = traitNames[match(H2s$Trait,traits)]
# H2s = subset(H2s,!is.na(Experimento) & TraitName %in% traitNames)
#
#
# files = list.files(path = 'single_trait_env_GP',pattern = 'GP_results',full.names=F)
#
# relative_results = do.call(bind_rows,lapply(files,function(file) {
#   trait = strsplit(file,'_')[[1]][3]
#   transform = sub('.csv','',strsplit(file,'_')[[1]][4])
#   res = read.csv(file.path('single_trait_env_GP',file),row.names=1)
#   res = do.call(bind_rows,tapply(1:nrow(res),res$Experimento,function(x)
#     data.frame(Trait = trait,Transform = transform,Experimento = res$Experimento[x[1]],PC_K = mean(res$PC_K[x]),t(colMeans(res[x,-c(1:2)]-res$PC_K[x])))))
# }))
#
# relative_results$Transform[relative_results$Transform == 'IND'] = 'INT'
# relative_results$h2 = H2s$h2[match(paste(relative_results$Experimento,relative_results$Trait),paste(H2s$Experimento,H2s$TraitName))]
# relative_results = relative_results[,!colnames(relative_results) == 'PC_K.1']
# relative_results_tall = tidyr::pivot_longer(relative_results,cols = -c(1:3,match(c('h2','PC_K'),colnames(relative_results))))
#
# results = do.call(bind_rows,lapply(files[1:10],function(file) {
#   trait = strsplit(file,'_')[[1]][3]
#   transform = sub('.csv','',strsplit(file,'_')[[1]][4])
#   res = read.csv(file.path('single_trait_env_GP',file),row.names=1)
#   res = do.call(bind_rows,tapply(1:nrow(res),res$Experimento,function(x)
#     data.frame(Trait = trait,Transform = transform,Experimento = res$Experimento[x[1]],t(colMeans(res[x,-c(1:2)])))))
# }))
#
#
# results$h2 = H2s$h2[match(paste(results$Experimento,results$Trait),paste(H2s$Experimento,H2s$TraitName))]
# results$h2[results$h2 < 0] = 0
#
# results$Transform[results$Transform == 'IND'] = 'INT'
# results_tall = tidyr::pivot_longer(results[,colnames(results) != 'h2'],cols = -c(1:3))
#
#
# library(ggplot2)
#
# pdf('GP_results.pdf')
# ggplot(results,aes(x=h2,y=PC_K)) + geom_point() + facet_grid(Transform~Trait)
# ggplot(results_tall,aes(x=name,y=value)) + geom_boxplot(aes(group = name,color = name)) + facet_grid(Trait~Transform,scales = 'free_y')
# ggplot(relative_results_tall,aes(x=name,y=log2((value+PC_K)/PC_K))) + geom_boxplot(aes(group = name,color = name)) + facet_grid(Trait~Transform,scales = 'free_y')
# dev.off()


results_long = pivot_longer(results, cols = c(tester, PC5, PC10, env, envPC, lm_K, lm_PC5, lm_env), names_to = 'model')
ggplot(results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'Experimento') +
  ylab('Predictive accuracy (r^2) improvement') +
  ggtitle('5-fold models improvement over tester variable model')

tester_results_long = pivot_longer(tester_results, cols = c(rf_tester, rf_PC5, rf_PC10, rf_env, rf_envPC, lm_K, lm_PC5, lm_env), names_to = 'model')
ggplot(tester_results_long, aes(x = model, y = value, color = model)) +
  geom_violin() +
  facet_wrap(facets = 'Experimento') +
  ylab('Predictive accuracy (r^2)') +
  ggtitle('Tester fold models improvement over tester variable model')
