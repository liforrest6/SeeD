library(MegaLMM)
library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(foreach)
library(doParallel)
registerDoParallel(RcppParallel::defaultNumThreads()-1)

run = commandArgs(t=T)
# run = as.numeric(strsplit(as.character(run),'')[[1]])
rep = as.numeric(run[1])
dataset = as.character(run[2])
dir = '/group/runciegrp2/Projects/SeeD/'

# rep = 1
# dataset = 'clim'
# traitN = 2

# print(c(floor(run/100),floor((run/10) %% 10),run %% 10))
print(c(rep,dataset))
# }
if(dataset == 'clim') {
  data = fread(paste0(dir, 'Env_data/GEA-climate-invnormtransformed.csv'), data.table = F)
# } else if (dataset == 'clim_log') {
#   data = fread('../data/SEEDGWAS_passport_environ_19clim_pheno_normalized.csv',data.table = F)
# } else if (dataset == 'PC_log') {
#   data = fread('../data/SEEDGWAS_passport_PCfactors_normalized.csv',data.table = F)
} else {
  data = NULL
}



# trait = colnames(data)[traitN]
results_folder = paste0(dir, sprintf('/Analyses/MegaLMM_output/dataset_%s-rep_%02d',dataset,rep))
try(dir.create(results_folder,recursive=T))

print(results_folder)

# data_wide = pivot_wider(data,id_cols = c('SampleID','Tester'),names_from = 'Experimento',values_from = all_of(trait))
# data_wide = subset(data_wide,rowSums(!is.na(data_wide[,-c(1:2)]))>0)
# data_wide = data_wide[,colSums(!is.na(data_wide))>0]
# 
# Y = as.matrix(data_wide[,-c(1:2)])
# Y = Y[,apply(Y,2,sd,na.rm=T)>1e-4]
# if(dataset == 3) {
#   Y = scale(Y)
# }
# 
# data_wide_std = cbind(data_wide[,1:2],Y)
# data_tall = pivot_longer(data_wide_std,cols = -c(1:2),names_to = 'Experimento',values_to = 'y')
# data_tall = subset(data_tall,!is.na(y))
# fwrite(data_tall,file = sprintf('%s/prepped_data.csv',results_folder))

## full data
Y = data[, -c(1)]
## test subsample for 100
# Y = data[1:100,-c(1)]
# Y = scale(Y)

sample_names = as.data.frame(data[, 1])
# subsample for 100
# sample_names = as.data.frame(data[1:100, 1])

colnames(sample_names) = c('SampleID')

K = fread('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

# subsample for test
# K = K[sample_names$SampleID,sample_names$SampleID]

k = ncol(Y)+1
run_parameters = MegaLMM_control(
  drop0_tol = 1e-10,#diagonalize_ZtZ_Kinv=F,
  scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 00,  # number of burn in samples before saving posterior samples
  K = k # number of factors
)


priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 18/20, nu = 100000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

runID = paste0(dir, sprintf('/Analyses/MegaLMM_output/MegaLMM_runs/%s_%d',dataset,rep))
MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                    ~ (1|SampleID),
                                    data=sample_names,         # the data.frame with information for constructing the model matrices
                                    relmat = list(SampleID = K),
                                    run_parameters=run_parameters,
                                    run_ID = runID
)
#maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y),verbose=T)
#map = maps$Missing_data_map_list[[which(maps$map_results$total_kept_NAs < 0.01*length(Y))[1]]]
#MegaLMM_state = set_Missing_data_map(MegaLMM_state,map) # keep < 1% of NAs
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
saveRDS(MegaLMM_state,file = sprintf('%s/MegaLMM_state_base.rds',runID))
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state = initialize_MegaLMM(MegaLMM_state)

# MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')
MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec')#,'Eta_mean')
# MegaLMM_state$Posterior$posteriorSample_params = c('Lambda')
MegaLMM_state$Posterior$posteriorMean_params = c()
MegaLMM_state$Posterior$posteriorFunctions = list(
  Eta_mean = 'F %*% Lambda + Z %*% U_R',
  G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
  R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])')
MegaLMM_state = clear_Posterior(MegaLMM_state)

n_iter = 200;  # how many samples to collect at once?
for(i  in 1:20) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  # print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
  
  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 10) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    print(MegaLMM_state$run_parameters$burn)
  }
}

G_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"G"))
R_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"R"))
P_hat = cov(Y,use='p')
write.csv(G_hat,file = sprintf('%s/G_hat.csv',results_folder))
write.csv(R_hat,file = sprintf('%s/R_hat.csv',results_folder))
write.csv(P_hat,file = sprintf('%s/P_hat.csv',results_folder))

