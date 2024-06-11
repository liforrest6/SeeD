library(data.table)
library(JointGWAS)
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)

run = commandArgs(t=T)
dataset = as.character(run[1])
chr_start = as.integer(run[2])
chr_end = as.integer(run[3])

dir = '/group/runciegrp2/Projects/SeeD/'

# set ncores. I use the RcppParallel function below, but you can set manually.
ncores = RcppParallel::defaultNumThreads()-1
registerDoParallel(ncores)

# old_genotype_dir = '/group/runciegrp2/Projects/SeeD_GWAS_Gates/Genetic_data/Beagle5/Genotype_Probabilities'
genotype_dir = paste0(dir, 'Genetic_data/Imputed_V4/genotypes_by_chromosome/')
if (dataset == 'clim') {
  output_dir = paste0(dir, 'Analyses/GEA_output/multivariate_results/') # set folder for output
} else if (dataset == 'clim_uni') {
  output_dir = paste0(dir, 'Analyses/GEA_output/univariate_results/') # set folder for output
}

print(output_dir)
# load environment data into variable Y
# should be the same as the input to MegaLMM
# should be no missing values
# rownames should be genotypeIDs
# column names should be environmental variables
# data = fread('as',data.table = F)

if (dataset == 'clim' | dataset == 'clim_perm') {
  data = as.data.frame(fread(file.path(dir, 'Env_data/GEA-climate-invnormtransformed.csv'),data.table = F))
  # load output from MegaLMM
# need matrices G_hat and R_hat
  G_hat = as.matrix(read.csv(file.path(dir, 'Analyses/MegaLMM_output/dataset_clim-rep_02/G_hat.csv'), row.names = 1))
  R_hat = as.matrix(read.csv(file.path(dir, 'Analyses/MegaLMM_output/dataset_clim-rep_02/R_hat.csv'), row.names = 1))

}


## for subsetting Y
# data = data[1:100, ]
# rownames(data) = data[1:100, 1]

# if (dataset == 'clim' | dataset == 'clim_perm') {
#   Y[, 2:8] = scale(Y[,2:8])
# }

SampleIDs = data[, 1]
# print(Y)

Y = data

dimnames(G_hat) = NULL
dimnames(R_hat) = NULL

# check that colnames(Y) == colnames(G_hat) == colnames(R_hat) == rownames(G_hat) == rownames(R_hat)

# run GWAS for each chromosome

# load genotype ID mapping file for GWAS columns
# TC_info = fread('/group/runciegrp2/Projects/SeeD_GWAS_Gates/prepped_data/TestCross_passport.csv',data.table=F)
non_dup_acc = read.csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/selected_genotypeIDs.csv')
colnames(non_dup_acc) = c('Unique.ID', 'SampleID')


for(chr in chr_start:chr_end) {
  print(chr)
  
  # load genotypes
  
  mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
  rownames(mat) = mat[,1]
  mat = as.matrix(mat[,-1])
  mat = mat[match(SampleIDs,non_dup_acc$Unique.ID),]
  rownames(mat) = SampleIDs
  # for permutation purposes, shuffle rownames of markers
  # rownames(mat) = sample(SampleIDs)


  print('finish loading genotypes')
  
  # discard genotypes with load maf or too much missing data (imputed)
  mean_notNA = colMeans(apply(mat,2,function(x) x %in% 0:2))
  maf = .5 - abs(.5-colMeans(mat)/2)
  select_cols = maf > 0.01 & mean_notNA>0.25
  mat = mat[,select_cols]

  print('finish filtering genotypes')
  
  # form map info for genotypes per chromosome (based on mat)
  map = fread(paste0(dir, 'Genetic_data/Imputed_V4/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.biallelic.012.pos'),
    data.table=F,
    col.names = c('V1', 'V2'))
  map$snp = sprintf('S%d_%d',map$V1,map$V2)
  map = map[match(colnames(mat),map$snp),]
  colnames(map) = c('Chr','pos','snp')
  map$maf = maf[select_cols]
  map$mean_imputed = 1-mean_notNA[select_cols]

  print('finish forming map')
  
  # load chromosome-specific K matrix
  K = fread(sprintf('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/LOO_Ks/K_chr_%02d.csv',chr),data.table=F)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  K = K[SampleIDs,SampleIDs]

  
  print('finish loading K matrix')

  # prep covariance matrices
  sK = svd(K)
  n = nrow(Y)
  sK2 = simultaneous_diagonalize(K,diag(1,n))
  sGR = simultaneous_diagonalize(G_hat,R_hat)

  print('finish calcing cov matrices')
  
  # As a test, to check that everything is working, 
  # subset to just the first 10 markers:
  # mat = mat[,1:100]
  # map = map[1:100,]


  results = EMMAX_ANOVA_matrix(cbind(tmin, tmax, trange, precipTot, aridityMean, rhMean, elevation)~X,
                            Y,mat,'Unique.ID',svd_matrices = list(sK,sGR),
                            mc.cores = ncores,verbose=T)  


  write.csv(cbind(map,results$anova),file = sprintf('%s/envGWAS_results_chr_%02d.csv',output_dir,chr),row.names=F)
  write.csv(results$beta_hats,file = sprintf('%s/envGWAS_beta_hats_chr_%02d.csv',output_dir,chr),row.names=F)
  write.csv(results$SEs,file = sprintf('%s/envGWAS_SEs_chr_%02d.csv',output_dir,chr),row.names=F)
}

print('Finished!')

