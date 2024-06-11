library(vcfR)
library(data.table)
library(tidyr)
library(dplyr)
library(JointGWAS)

library(Matrix)

genotype_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome'
LOO_K_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/LOO_Ks'
gemma = '~/software/gemma/gemma-0.98.5-linux-static-AMD64'

run = as.numeric(commandArgs(t=T)[1])
run = as.numeric(strsplit(as.character(run),'')[[1]])

rep = run[1]
dataset = run[2]
transform = run[3]
traitN = run[4]
chr = run[5]
if(chr == 0) chr = 10
if(is.na(rep)) rep=1
if(is.na(dataset)) dataset=1
if(is.na(transform)) transform=1
if(is.na(traitN)) traitN = 4
if(is.na(chr)) chr = 10
if(chr == 0) chr=10

print(c(rep,dataset,transform,traitN,chr))
trait = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")[traitN]
results_folder = sprintf('%s/dataset_%02d_%s/rep_%02d/GEMMA_univariate_byTester',trait,dataset,ifelse(transform==1,'orig','INT'),rep)
try(dir.create(results_folder,recursive=T))

# if(file.exists(sprintf('%s/GxE_GWAS_results_chr_%02d.csv',results_folder,chr))) q()

data_tall = fread(file.path(results_folder,'../prepped_data.csv'),data.table=F)
data_wide = pivot_wider(data_tall,id_cols = c('SampleID'),names_from = c('Experimento','Tester'),values_from = 'y')
colnames(data_wide) = gsub('/',':',colnames(data_wide))
Y = as.matrix(data_wide[,-c(1)])
rownames(Y) = data_wide$SampleID

Y = Y[,apply(Y,2,sd,na.rm=T)>1e-4 & colSums(!is.na(Y)) >= 50]
if(dataset == 3) {
  Y = scale(Y)
}
if(transform == 2) {
  Y = apply(Y,2,function(x) {
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  })
}

# remove section when cholL_Sigma_inv already built
SampleIDs = data_wide$SampleID
mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
rownames(mat) = mat[,1]
mat = as.matrix(mat[,-1])
mat = mat[SampleIDs,]


map = fread(file.path(genotype_dir,'../V2_V4_mapping.csv'),data.table=F)
map = map[match(colnames(mat),map$V4),]

setwd(results_folder)


K_file = sprintf('%s/cholL_Sigma_inv/cholL_Sigma_inv_chr%02d.txt',results_folder,chr)
K = fread(sprintf('%s/K_chr_%02d.csv',LOO_K_dir,chr),data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
K = K[SampleIDs,SampleIDs]

MAF = beta_hat = SE = p_wald = matrix(NA,nrow = nrow(map),ncol = ncol(Y),dimnames = list(map$V4,colnames(Y)))

for(experimento in colnames(Y)) {
  print(sprintf('column %d of %d',match(experimento,colnames(Y)),ncol(Y)))
  y = Y[,experimento]
  i = !is.na(y)
  X = mat[i,]
  y = y[i]
  Ki = K[i,i]

  maf = .5 - abs(.5-colMeans(X)/2)
  select_cols = maf > 0.01

  X = X[,select_cols]
  map_exp = map[match(colnames(X),map$V4),]

  geno = data.frame(ID = colnames(X),REF='A',ALT='T',t(X))

  tmp_file = paste(sample(letters,5),collapse='')

  fwrite(geno,file = sprintf('%s_geno.bimbam',tmp_file),row.names=F,col.names = F,sep=',')
  write.table(y,file = sprintf('%s_pheno.txt',tmp_file),sep='\t',row.names=F,col.names=F)
  write.table(Ki,file = sprintf('%s_K.txt',tmp_file),row.names=F,col.names=F)

  outfile = sprintf('gemma_output_%s_chr%02d.txt',experimento,chr)

  system(sprintf('%s -g %s -p %s -k %s -lmm 4 -o %s',
                 gemma,sprintf('%s_geno.bimbam',tmp_file),sprintf('%s_pheno.txt',tmp_file),sprintf('%s_K.txt',tmp_file),
                 outfile))
  system(sprintf('rm -rf %s*',tmp_file))
  results = fread(sprintf('output/%s.assoc.txt',outfile),data.table=F)
  MAF[select_cols,experimento] = results$af
  beta_hat[select_cols,experimento] = results$beta
  SE[select_cols,experimento] = results$se
  p_wald[select_cols,experimento] = results$p_wald
}

i = rowSums(!is.na(MAF))>0

fwrite(as.data.frame(MAF[i,]),file = sprintf('MAF_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(beta_hat[i,]),file = sprintf('beta_hat_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(SE[i,]),file = sprintf('SE_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')
fwrite(as.data.frame(p_wald[i,]),file = sprintf('p_wald_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')


