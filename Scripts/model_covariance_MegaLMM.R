library(MegaLMM)
library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(foreach)
library(doParallel)
registerDoParallel(RcppParallel::defaultNumThreads()-1)
MegaLMM::set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)

run = as.numeric(commandArgs(t=T)[1])
run = as.numeric(strsplit(as.character(run),'')[[1]])
rep = run[1]
dataset = run[2]
transform = run[3]
traitN = run[4]

# print(c(floor(run/100),floor((run/10) %% 10),run %% 10))
print(c(rep,dataset,traitN))
# }
if(dataset == 1) {
  data = fread('../../Phenotype_data/blups_deregressed.csv',data.table = F)
} else if(dataset == 2) {
  data = fread('../../Phenotype_data/blups_deregressed_CV.csv',data.table = F)
} else if(dataset == 3) {
  data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
}

data = subset(data,Trial_Classification != 'Stress')

trait = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")[traitN]

results_folder = sprintf('%s/dataset_%02d_%s/rep_%02d',trait,dataset,ifelse(transform==1,'orig','INT'),rep)
try(dir.create(results_folder,recursive=T))

data_wide = pivot_wider(subset(data,Trait == trait),id_cols = c('SampleID','Tester'),names_from = 'Experimento',values_from = 'Value')
data_wide = subset(data_wide,rowSums(!is.na(data_wide[,-c(1:2)]))>0)
data_wide = data_wide[,colSums(!is.na(data_wide))>0]

sample_to_geneticData = fread('../../Genetic_data/Imputed_V4/selected_genotypeIDs.csv',data.table=F)
data_wide = subset(data_wide,SampleID %in% sample_to_geneticData$Sample)
data_wide$SampleID = sample_to_geneticData$V1[match(data_wide$SampleID,sample_to_geneticData$Sample)]

Y = as.matrix(data_wide[,-c(1:2)])
Y = Y[,apply(Y,2,sd,na.rm=T)>1e-4]
if(transform == 2) {
  Y = apply(Y,2,function(x) {
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  })
}


Y_resid = sapply(1:ncol(Y),function(i) {
  j = !is.na(Y[,i])
  r = resid(lm(Y[,i]~Tester,data_wide))
  y = rep(NA,nrow(Y))
  y[j] = r
  y})
colnames(Y_resid) = colnames(Y)
Y = Y_resid

# Y = Y[,c(3,21)]
# j = !is.na(rowSums(Y))
# Y = Y[j,]
# # Y = cbind(Y[j,3],Y[j,3])
# data_wide = data_wide[j,]

data_wide_std = cbind(data_wide[,1:2],Y)
data_tall = pivot_longer(data_wide_std,cols = -c(1:2),names_to = 'Experimento',values_to = 'y')
data_tall = subset(data_tall,!is.na(y))
fwrite(data_tall,file = sprintf('%s/prepped_data.csv',results_folder))

K = fread('../../Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
K = K[data_wide$SampleID,data_wide$SampleID]

# make K for FOAM
diag(K) = 1
K = K/4

#
# Y = Y[1:1000,]
# data_wide = data_wide[1:1000,]
# K = K[data_wide$SampleID,data_wide$SampleID]
# Y = Y[,colSums(!is.na(Y))>100]
#
# var_comp = do.call(cbind,mclapply(1:ncol(Y),function(i) {
#   j = !is.na(Y[,i])
#   X = model.matrix(~Tester,data_wide[j,])
#   res = mixed.solve(Y[j,i],K=K[j,j],X=X)
#   c(Vu = res$Vu,Ve = res$Ve)
# },mc.cores = RcppParallel::defaultNumThreads()))
# #
# us = do.call(cbind,mclapply(1:ncol(Y),function(i) {
#   j = !is.na(Y[,i])
#   X = model.matrix(~Tester,data_wide[j,])
#   res = mixed.solve(Y[j,i],K=K[j,j],X=X)
#   u = rep(NA,nrow(Y))
#   u[j] = res$u
#   u
# },mc.cores = RcppParallel::defaultNumThreads()))

k = ncol(Y)-1
n_h2 = 30
run_parameters = MegaLMM_control(
  drop0_tol = 1e-10,#diagonalize_ZtZ_Kinv=F,
  scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = n_h2, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 00,  # number of burn in samples before saving posterior samples
  K = k # number of factors
)


priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.1,   nu = 2),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 1, nu = 100000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_horseshoe,
  #   prop_0 = 0.1,
  #   delta = list(shape = 3, scale = 1),
  #   delta_iterations_factor = 100
  # ),
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    prop_0 = 0.1,
    Lambda_df   = 2,
    delta_1 = list(shape = 2.1,rate = 1/2),
    delta_2 = list(shape = 3,rate = 2),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) ifelse(h2s==(n-1)/n,n,n/(n-1)) # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

runID = sprintf('MegaLMM_runs/%s_%d_%d',trait,dataset,rep)
MegaLMM_state = setup_model_MegaLMM(Y_resid,            # n x p data matrix
                                    ~ Tester + (1|SampleID),
                                    data=data_wide,         # the data.frame with information for constructing the model matrices
                                    relmat = list(SampleID = K),
                                    run_parameters=run_parameters,
                                    run_ID = runID
)
maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y)+1,verbose=T)
map = maps$Missing_data_map_list[[which(maps$map_results$total_kept_NAs < 0.01*length(Y))[1]]]
MegaLMM_state = set_Missing_data_map(MegaLMM_state,map) # keep < 1% of NAs
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
saveRDS(MegaLMM_state,file = sprintf('%s/MegaLMM_state_base.rds',runID))
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state = initialize_MegaLMM(MegaLMM_state)

# MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')
MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec','B1')#,'Eta_mean')
# MegaLMM_state$Posterior$posteriorSample_params = c('Lambda')
MegaLMM_state$Posterior$posteriorFunctions = list(
  G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
  R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
  h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])')
MegaLMM_state = clear_Posterior(MegaLMM_state)

n_iter = 200;  # how many samples to collect at once?
for(i  in 1:20) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 10) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    print(MegaLMM_state$run_parameters$burn)
  }
}
# data_wide$y = Y[,1]
# res = GridLMM::GridLMM_posterior(y~Tester+(1|SampleID),data = data_wide,relmat = list(SampleID = K),h2_divisions = 200,normalize_relmat = F)
# post = res$h2s_results
# post = post[order(post[,1]),]
h2s = load_posterior_param(MegaLMM_state,"h2")
saveRDS(h2s,file = sprintf('%s/h2s.csv',results_folder))
# plot(post[,1],post[,2]/diff(post[,1])[1],type='l')
# lines(density(h2s[,1,1]),col=2)
# boxplot(load_posterior_param(MegaLMM_state,'F_h2')[,1,])
# boxplot(load_posterior_param(MegaLMM_state,'resid_h2')[,1,])
# MegaLMM_state$Posterior$Lambda = load_posterior_param(MegaLMM_state,'Lambda')
# MegaLMM_state$Posterior$tot_Eta_prec = load_posterior_param(MegaLMM_state,'tot_Eta_prec')
# boxplot(get_posterior_FUN(MegaLMM_state,colSums(Lambda^2)*tot_Eta_prec[1,]))
# hist(h2s[,1,1]);hist(h2s[,2,1]);dev.off()
G_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"G"))
R_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"R"))
# cov2cor(G_hat)
# R_hat
# colMeans(h2s[,,1])
# diag(G_hat)/diag(G_hat+R_hat)

write.csv(G_hat,file = sprintf('%s/G_hat.csv',results_folder))
write.csv(R_hat,file = sprintf('%s/R_hat.csv',results_folder))

# Y_resid = sapply(1:ncol(Y),function(i) {
#   j = !is.na(Y[,i])
#   r = resid(lm(Y[,i]~Tester,data_wide))
#   y = rep(NA,nrow(Y))
#   y[j] = r
#   y})
# colnames(Y_resid) = colnames(Y)
P_hat = cov(Y_resid,use='p')
write.csv(P_hat,file = sprintf('%s/P_hat.csv',results_folder))

# G_hat1 = get_posterior_mean(load_posterior_param(MegaLMM_state,"G")[1:100,,])
# G_hat2 = get_posterior_mean(load_posterior_param(MegaLMM_state,"G")[900+1:100,,])
