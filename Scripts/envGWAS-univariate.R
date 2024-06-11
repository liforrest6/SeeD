library(vcfR)
library(data.table)
library(tidyr)
library(dplyr)

# library(JointGWAS)

library(Matrix)

genotype_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome'
LOO_K_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/LOO_Ks'
gemma = '~/software/gemma/gemma-0.98.5-linux-static-AMD64'
dir = '/group/runciegrp2/Projects/SeeD/'

run = as.numeric(commandArgs(t=T)[1])
run = as.numeric(strsplit(as.character(run),'')[[1]])

traitN = run[1]
rep = run[2]
dataset = run[3]
transform = run[4]

# chr = run[5]
if(is.na(traitN)) traitN = 4
if(is.na(rep)) rep=1
if(is.na(dataset)) dataset=1
if(is.na(transform)) transform=1
# if(is.na(chr)) chr = 4
# if(chr == 0) chr = 4


# print(c(rep,dataset,transform,traitN,chr))
print(c(rep,dataset,transform,traitN))
trait = c('tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')[traitN]
print(trait)
print(paste('rep:', rep))
results_folder = sprintf('%s/Analyses/GEA_output/GEMMA_univariate/%s/rep_%02d/',dir,trait,rep)
try(dir.create(results_folder,recursive=T))

# if(file.exists(sprintf('%s/GxE_GWAS_results_chr_%02d.csv',results_folder,chr))) q()

data = fread(file.path(dir,'Env_data', 'GEA-climate-invnormtransformed.csv'),data.table=F)
# data_wide = pivot_wider(data_tall,id_cols = c('SampleID','Tester'),names_from = 'Experimento',values_from = 'y')
# Y = as.matrix(data_wide[,-c(1:2)])

# Y = Y[,apply(Y,2,sd,na.rm=T)>1e-4]
# if(dataset == 3) {
#   Y = scale(Y)
# }
# if(transform == 2) {
#   Y = apply(Y,2,function(x) {
#     qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
#   })
# }

# testing with 20 accessions
# data = data[1:20,]
Y = data

for(chr in 1:10){

  # remove section when cholL_Sigma_inv already built
  SampleIDs = data$Unique.ID
  mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
  rownames(mat) = mat[,1]
  mat = as.matrix(mat[,-1])
  mat = mat[SampleIDs,]

  # testing with 200 SNPs
  # mat = mat[SampleIDs, 1:200]


  map = fread('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/V2_V4_mapping.csv',data.table=F)
  map = map[match(colnames(mat),map$V4),]

  setwd(results_folder)


  # K_file = sprintf('%s/cholL_Sigma_inv/cholL_Sigma_inv_chr%02d.txt',results_folder,chr)
  K = fread(sprintf('%s/K_chr_%02d.csv',LOO_K_dir,chr),data.table=F)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  K = K[SampleIDs,SampleIDs]

  MAF = beta_hat = SE = p_wald = matrix(NA,nrow = nrow(map),ncol = ncol(Y),dimnames = list(map$V4,colnames(Y)))

# for(trait in colnames(Y)) {
  print(sprintf('column %d of %d',match(trait,colnames(Y)),ncol(Y)))
  y = Y[,trait]
  i = !is.na(y)
  # print(y)
  # print(length(y))
  # print(length(i))
  # print(mat)
  # print(length(mat))

  X = mat[i,]
  y = y[i]
  Ki = K[i,i]
  # Xt = model.matrix(~,data[i,])

  maf = .5 - abs(.5-colMeans(X)/2)
  select_cols = maf > 0.01

  X = X[,select_cols]
  map_exp = map[match(colnames(X),map$V4),]

  geno = data.frame(ID = colnames(X),REF='A',ALT='T',t(X))

  tmp_file = paste(sample(letters,5),collapse='')
  print(tmp_file)

  fwrite(geno,file = sprintf('%s_geno.bimbam',tmp_file),row.names=F,col.names = F,sep=',')
  # write.table(Xt,file = sprintf('%s_cov.txt',tmp_file),row.names=F,col.names = F,sep=',')
  write.table(y,file = sprintf('%s_pheno.txt',tmp_file),sep='\t',row.names=F,col.names=F)
  write.table(Ki,file = sprintf('%s_K.txt',tmp_file),row.names=F,col.names=F)

  outfile = sprintf('gemma_output_chr%02d.txt',chr)

  # system(sprintf('%s -g %s -p %s -k %s  -c %s -lmm 4 -o %s',
  #                'gemma',sprintf('%s_geno.bimbam',tmp_file),sprintf('%s_pheno.txt',tmp_file),sprintf('%s_K.txt',tmp_file),sprintf('%s_cov.txt',tmp_file),
  #                outfile))
  ## testing without cov matrix, with no tester fixed efect
  system(sprintf('%s -g %s -p %s -k %s  -lmm 4 -o %s',
               'gemma',sprintf('%s_geno.bimbam',tmp_file),sprintf('%s_pheno.txt',tmp_file),sprintf('%s_K.txt',tmp_file),
               outfile))
  system(sprintf('rm -rf %s*',tmp_file))
  results = fread(sprintf('output/%s.assoc.txt',outfile),data.table=F)
  MAF[select_cols,trait] = results$af
  beta_hat[select_cols,trait] = results$beta
  SE[select_cols,trait] = results$se
  p_wald[select_cols,trait] = results$p_wald
# }

i = rowSums(!is.na(MAF))>0

fwrite(as.data.frame(MAF[i,]),file = sprintf('MAF_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(beta_hat[i,]),file = sprintf('beta_hat_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(SE[i,]),file = sprintf('SE_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(p_wald[i,]),file = sprintf('p_wald_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')

}
