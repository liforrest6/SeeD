library(MegaLMM)
library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(spatialsample)
library(sf)
library(SpatialKDE)
library(sp)
library(tmap)
library(foreach)
library(doParallel)
library(caret)
registerDoParallel()

run = commandArgs(t=T)
# run = as.numeric(strsplit(as.character(run),'')[[1]])
rep = as.numeric(run[1])
i_fold = as.numeric(run[2])

dir = '/group/runciegrp2/Projects/SeeD/'

# rep = 2
dataset = 'clim'
# traitN = 2

# print(c(floor(run/100),floor((run/10) %% 10),run %% 10))
print(c(rep, i_fold))

clusters = read.csv(file.path(dir, sprintf('Analyses/MegaLMM_output/dataset_clim-spatial_prediction-rep_%02d/testing_folds_assignment.csv', rep)))
fold_accessions = clusters %>% filter(fold == sprintf('Fold%d', i_fold))

K = fread('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

env_data = fread(paste0(dir, 'Env_data/GEA-climate-invnormtransformed.csv'), data.table = F)
geno_info = read.csv(file.path(dir, 'Phenotype_data/selected_genotypeIDs.csv'))
passport_data = readxl::read_excel(file.path(dir, 'Env_data/Passport/Maize climate data nov 2022 AEZ elevation by growing season and flowering period.xlsx'))
dna_to_tc_gid = as.data.frame(readxl::read_xlsx(file.path(dir, 'Phenotype_data/Hearne_data/dataset/Blups_01.xlsx'),sheet = 'DNA to TC GID'))


geno_info$AccID = dna_to_tc_gid[match(geno_info$Sample,dna_to_tc_gid$`Sample ID`),2]
geno_info$LatNew = passport_data$LatNew[match(geno_info$AccID,passport_data$GID)]
geno_info$LongNew = passport_data$LongNew[match(geno_info$AccID,passport_data$GID)]

geno_info = geno_info[!is.na(geno_info$LatNew) & geno_info$V1 %in% colnames(K) & geno_info$V1 %in% env_data$Unique.ID,]
K = K[geno_info$V1,geno_info$V1]
env_data = env_data[match(geno_info$V1,env_data$Unique.ID),]

geno_info = merge(geno_info,env_data,by.x = 'V1',by.y = 'Unique.ID')

geno_sf = st_as_sf(geno_info,coords = c('LongNew','LatNew'))
# need to set coordinate system
st_crs(geno_sf) <- 4326

## full data

## test subsample for 100
# Y = data[1:100,-c(1)]
# Y = scale(Y)

# trait = colnames(data)[traitN]
results_folder = paste0(dir, sprintf('/Analyses/MegaLMM_output/dataset_%s-spatial_prediction-rep_%02d',dataset,rep))
# try(dir.create(results_folder,recursive=T))

print(results_folder)

# combined_results = foreach(i_fold = 1:nrow(clusters),.combine=rbind) %do%   {

  print(sprintf('Start fold %d', i_fold))

  fold_data = env_data

  testing_index = match(fold_accessions$accession,fold_data$Unique.ID)
  # training_index = match(training$V1,fold_data$Unique.ID)

  
  validationIDs = env_data[testing_index,]$Unique.ID
  fold_data[testing_index, -c(1)] = NA

  # sampled_data = data[folds[[i_fold]],]
  # validationIDs = sampled_data$Unique.ID
  # Y[which(data$Unique.ID %in% validationIDs),] = NA
  write.csv(as.data.frame(validationIDs), file = sprintf('%s/withheld-accessions-rep_%02d-fold_%02d.csv',results_folder, rep, i_fold))


  sample_names = as.data.frame(env_data[, 1])
  # subsample for 100
  # sample_names = as.data.frame(data[1:100, 1])

  colnames(sample_names) = c('SampleID')
  Y_mat = fold_data[, -c(1)]

  k_num = ncol(Y_mat)+1


  run_parameters = MegaLMM_control(
    drop0_tol = 1e-10,#diagonalize_ZtZ_Kinv=F,
    scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 00,  # number of burn in samples before saving posterior samples
    K = k_num # number of factors
  )

  print('Finish parameter set-up')


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

  runID = paste0(dir, sprintf('/Analyses/MegaLMM_output/MegaLMM_runs/%s-spatial_validation_%d_%d',dataset,rep, i_fold))
  MegaLMM_state = setup_model_MegaLMM(Y_mat,            # n x p data matrix WITH INFORMATION WITHHELD
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
  saveRDS(MegaLMM_state,file = sprintf('%s/MegaLMM_state_base_%d.rds',runID, i_fold))
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state)

  print('Finish MegaLMM initialization')

  # MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')
  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec', 'U_F', 'U_R')#,'Eta_mean')
  # MegaLMM_state$Posterior$posteriorSample_params = c('Lambda')
  MegaLMM_state$Posterior$posteriorMean_params = c()
  MegaLMM_state$Posterior$posteriorFunctions = list(
    Eta_mean = 'F %*% Lambda + Z %*% U_R',
    G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
    R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])')
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  print(MegaLMM_state$Posterior$folder)
  print(MegaLMM_state$Posterior$files)

  print('Begin MCMC')

  # n_iter = 200;  # how many samples to collect at once?
  n_iter = 200
  for(i  in 1:20) {
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
    
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
    
    print(MegaLMM_state$Posterior$folder)
    print(MegaLMM_state$Posterior$files)

    # set of commands to run during burn-in period to help chain converge
    if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 10) {
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
      print(MegaLMM_state$run_parameters$burn)
    }
  }
  print('Finish MCMC')

    print(MegaLMM_state$Posterior$folder)
  print(MegaLMM_state$Posterior$files)


  # G_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"G"))
  # R_hat = get_posterior_mean(load_posterior_param(MegaLMM_state,"R"))
  # P_hat = cov(Y,use='p')
  # write.csv(G_hat,file = sprintf('%s/G_hat-rep_%02d.csv',results_folder, rep))
  # write.csv(R_hat,file = sprintf('%s/R_hat-rep_%02d.csv',results_folder, rep))
  # write.csv(P_hat,file = sprintf('%s/P_hat-rep_%02d.csv',results_folder, rep))

  MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
  U_hat = get_posterior_mean(MegaLMM_state,'U_R + U_F %*% Lambda')
  write.csv(U_hat,file = sprintf('%s/U_hat-rep_%02d-fold_%02d.csv',results_folder, rep, i_fold))

  # U_hat_withheld = U_hat[testing_index,]
  print(sprintf('Finished with fold %d', i_fold))
  # U_hat_withheld
# }

# write.csv(combined_results,file = sprintf('%s/U_hat-rep_%02d-combined.csv',results_folder, rep))


## END

