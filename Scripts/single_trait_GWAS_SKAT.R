library(data.table)
# library(dplyr)
library(JointGWAS)

library(Matrix)
library(rrBLUP)

library(SKAT)
library(foreach)
library(doParallel)
registerDoParallel(4)

genotype_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome'
LOO_K_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/LOO_Ks'
gemma = '~/software/gemma/gemma-0.98.5-linux-static-AMD64'
gff_file= '/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.36.gff3'

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
results_folder = sprintf('%s/dataset_%02d_%s/rep_%02d',trait,dataset,ifelse(transform==1,'orig','INT'),rep)
try(dir.create(results_folder,recursive=T))

if(file.exists(sprintf('%s/SKAT_unweighted_p_chr%02d.txt',results_folder,chr))) q()

K = fread('../../Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
sample_to_geneticData = fread('../../Genetic_data/Imputed_V4/selected_genotypeIDs.csv',data.table=F)

data = fread('../../Phenotype_data/blups_std.csv',data.table = F)
trait = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")[traitN]
data$DNAID = sample_to_geneticData$V1[match(data$SampleID,sample_to_geneticData$Sample)]

data = subset(data,Trait == trait & !is.na(DNAID))

# remove section when cholL_Sigma_inv already built
SampleIDs = unique(data$DNAID)
mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
rownames(mat) = mat[,1]
mat = as.matrix(mat[,-1])
mat = mat[SampleIDs,]

map = fread(file.path(genotype_dir,'../V2_V4_mapping.csv'),data.table=F)
map = map[match(colnames(mat),map$V4),]

gff = fread(cmd=sprintf('grep "ID=gene:" %s',gff_file),data.table=F)
gff = subset(gff,V1==chr)
gff$Gene = sapply(gff$V9,function(x) sub('ID=gene:','',strsplit(x,';')[[1]][1]))
gff = gff[order(gff$V4),]
gff$cut = c((gff$V5[-nrow(gff)]+gff$V4[-1])/2,NA)

gene_sets = list()
i = 1
end = 0
while(i <= nrow(gff)) {
  if(i %% 1000 == 0) print(i)
  j = i
  start = end + 1
  end = gff$V5[j]
  if(j < nrow(gff)) {
    while(end > gff$V4[j+1]) {
      j = j+1
      if(j >= (nrow(gff)-2)) break
      end = max(end,gff$V5[j])
    }
    # if(j > i+1) break
    end = (end + gff$V4[j+1])/2
  }
  gene_sets[[length(gene_sets)+1]] = data.frame(Genes = paste(gff$Gene[i:j],collapse=','),start = start,end = end)
  i = j+1
  if(i > (nrow(gff))) break
}
gene_sets = do.call(rbind,gene_sets)
gene_sets$end[nrow(gene_sets)] = max(map$POS)

map$Gene = NA

gene_set_index = 1
map$Gene = NA
for(i in 1:nrow(map)) {
  while(map[i,2]>gene_sets$end[gene_set_index]) {
    gene_set_index = gene_set_index + 1
    if(gene_set_index > nrow(gene_sets)) {
      gene_set_index = nrow(gene_sets)
      break
    }
  }
  map$Gene[i] = gene_sets$Gene[gene_set_index]
}

# setwd(results_folder)


K_file = sprintf('%s/cholL_Sigma_inv/cholL_Sigma_inv_chr%02d.txt',results_folder,chr)
K = fread(sprintf('%s/K_chr_%02d.csv',LOO_K_dir,chr),data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
K = K[SampleIDs,SampleIDs]

results = foreach(exp = unique(data$Experimento),.combine = cbind) %dopar% {
  data_env = subset(data,Experimento == exp)
  if(transform == 2) {
    x = data_env$Value
    data_env$Value = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  }

  K_env = K[data_env$DNAID,data_env$DNAID]

  X = model.matrix(~Tester,data_env)

  obj<-SKAT_NULL_emmaX(data_env$Value ~ X, K=K_env)

  p = rep(NA,nrow(gene_sets))
  for(i in 1:nrow(gene_sets)) {
    gene = gene_sets$Gene[i]
    j = map$Gene == gene
    if(sum(j) == 0) next
    Z = mat[data_env$DNAID,j,drop=FALSE]
    Z = Z[,apply(Z,2,var)>0,drop=FALSE]
    if(ncol(Z)<1) next
    Z = apply(Z,2,function(x) {if(mean(x)>1) x = 2-x;x})
    if(ncol(Z)<1) next
    Z = Z[,colMeans(Z)/2 > 0.01,drop=FALSE]
    if(ncol(Z)<1) next
    res = try(SKAT(Z, obj,weights = rep(1,ncol(Z))))
    if(is(res,'SKAT_OUT')) {
      p[i] = res$p.value
    } else{
      p[i] = -1
    }
  }
  p
}
colnames(results) = unique(data$Experimento)
results = cbind(data.frame(CHROM=chr,gene_sets),results)
setwd(results_folder)
fwrite(results,file = sprintf('SKAT_unweighted_p_chr%02d.txt',chr),row.names=T,col.names=T,sep=',',na = 'NA')


