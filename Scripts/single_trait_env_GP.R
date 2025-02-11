# how well does genomic prediction work within trials for each trait?
# does it work better than PCs?
library(data.table)
library(dplyr)
library(rrBLUP)
library(foreach)
library(doParallel)
registerDoParallel(4)


run = as.numeric(commandArgs(t=T)[1])
run = as.numeric(strsplit(as.character(run),'')[[1]])
# rep = run[1]
# dataset = run[2]
transform = run[1]
traitN = run[2]

if(is.na(transform)) transform = 1
if(is.na(traitN)) traitN = 4

K = fread('../../Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
sample_to_geneticData = fread('../../Genetic_data/Imputed_V4/selected_genotypeIDs.csv',data.table=F)

data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
trait = c("ASI","DaysToFlowering","FieldWeight", "BareCobWeight","PlantHeight", "GrainWeightPerHectareCorrected","GrainWeightPerHectare")[traitN]
data$DNAID = sample_to_geneticData$V1[match(data$SampleID,sample_to_geneticData$Sample)]

data = subset(data,Trait == trait & !is.na(DNAID))

# results = foreach(exp = unique(data$Experimento),.combine = bind_rows) %dopar% {
#   data_env = subset(data,Experimento == exp)
#   if(transform == 2) {
#     x = data_env$Value
#     data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
#   }
#
#   K_env = K[data_env$DNAID,data_env$DNAID]
#   sK = svd(K_env)
#
#   X = model.matrix(~Tester,data_env)
#   X_PC = cbind(X,sK$u[,1:10])
#
#   folds = caret::createFolds(data_env$Value,k=5)
#   foreach(i_fold = seq_along(folds),.combine=bind_rows) %do% {
#     y = data_env$Value
#     y[folds[[i_fold]]] = NA
#
#     res_X = mixed.solve(y,X = X)
#     res_PC5 = mixed.solve(y,X = cbind(X,sK$u[,1:5]))
#     res_PC10 = mixed.solve(y,X = cbind(X,sK$u[,1:10]))
#     res_PC_K = mixed.solve(y,X = X_PC,K=K_env)
#     res_X_K = mixed.solve(y,X = X,K=K_env)
#
#     data.frame(Experimento = exp, fold = i_fold,
#                X = cor(cbind(data_env$Value,X %*% res_X$beta)[folds[[i_fold]],])[2],
#                PC5 = cor(cbind(data_env$Value,cbind(X,sK$u[,1:5]) %*% res_PC5$beta)[folds[[i_fold]],])[2],
#                PC10 = cor(cbind(data_env$Value,cbind(X,sK$u[,1:10]) %*% res_PC10$beta)[folds[[i_fold]],])[2],
#                PC_K = cor(cbind(data_env$Value,c(X_PC %*% res_PC_K$beta)+res_PC_K$u)[folds[[i_fold]],])[2],
#                X_K = cor(cbind(data_env$Value,c(X %*% res_X_K$beta) + res_X_K$u)[folds[[i_fold]],])[2])
#   }
# }
#
# write.csv(results,file = sprintf('single_trait_env_GP/GP_results_%s_%s.csv',trait,ifelse(transform == 1,'orig','IND')))


results = foreach(exp = unique(data$Experimento),.combine = bind_rows,.errorhandling = 'remove') %dopar% {
  data_env = subset(data,Experimento == exp)
  if(transform == 2) {
    x = data_env$Value
    data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  }

  K_env = K[data_env$DNAID,data_env$DNAID]
  sK = svd(K_env)

  X = model.matrix(~Tester,data_env)
  X_PC = cbind(X,sK$u[,1:10])

  foreach(tester = unique(data_env$Tester),.combine=bind_rows,.errorhandling = 'remove') %do% {
    y = data_env$Value
    testing = data_env$Tester == tester
    y[testing] = NA

    dropcols = caret::findLinearCombos(X[!testing,])
    X_tester = X
    if(length(dropcols$remove) > 0) {
      X_tester = X[,-dropcols$remove]
    }

    res_X = mixed.solve(y,X = X_tester)
    res_PC5 = mixed.solve(y,X = cbind(X_tester,sK$u[,1:5]))
    res_PC10 = mixed.solve(y,X = cbind(X_tester,sK$u[,1:10]))
    res_PC_K = mixed.solve(y,X = cbind(X_tester,sK$u[,1:10]),K=K_env)
    res_X_K = mixed.solve(y,X = X_tester,K=K_env)

    data.frame(Experimento = exp, Tester = tester,n_train = sum(!is.na(y)),n_test = sum(is.na(y)),
               X = cor(cbind(data_env$Value,X_tester %*% res_X$beta)[testing,])[2],
               PC5 = cor(cbind(data_env$Value,cbind(X_tester,sK$u[,1:5]) %*% res_PC5$beta)[testing,])[2],
               PC10 = cor(cbind(data_env$Value,cbind(X_tester,sK$u[,1:10]) %*% res_PC10$beta)[testing,])[2],
               PC_K = cor(cbind(data_env$Value,c(cbind(X_tester,sK$u[,1:10]) %*% res_PC_K$beta)+res_PC_K$u)[testing,])[2],
               X_K = cor(cbind(data_env$Value,c(X_tester %*% res_X_K$beta) + res_X_K$u)[testing,])[2])
  }
}

write.csv(results,file = sprintf('single_trait_env_GP/GP_tester_results_%s_%s.csv',trait,ifelse(transform == 1,'orig','INT')))



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



