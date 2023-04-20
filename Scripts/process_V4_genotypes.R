library(data.table)

vcf_file = '/group/runciegrp/SharedResources/CIMMYT/Imputed_v4/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf'

vcf = fread(cmd=sprintf('tail -n +11 %s | cut -f 1-8',vcf_file),data.table=F)

v2_v4 = data.frame(CHROM = vcf[,1],POS = vcf[,2],ID = vcf$ID,V4 = sprintf('S%s_%s',vcf[,1],vcf[,2]))

genotypes = fread('zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.biallelic.012',data.table=F)
genotypes = genotypes[,-1]
indiv = fread('zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.biallelic.012.indv',data.table=F,h=F)
indiv$Sample = sapply(indiv[,1],function(x) strsplit(x,':')[[1]][1])

pos = fread('zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.biallelic.012.pos',data.table=F,h=F)
pos$name = sprintf('S%s_%s',pos[,1],pos[,2])

# filter for positions that remained after biallelic filetering
v2_v4 = v2_v4[match(pos$name,v2_v4$V4),]
write.csv(v2_v4,file = 'V2_V4_mapping.csv',row.names=F)

# just choose the first of each sample
selected_individuals = indiv[match(unique(indiv$Sample),indiv$Sample),]
write.csv(selected_individuals,file = 'selected_genotypeIDs.csv',row.names=F)

genotypes = genotypes[match(selected_individuals$V1,indiv$V1),]
colnames(genotypes) = pos$name

# split by chromosome
try(dir.create('genotypes_by_chromosome'))
for(chr in 1:10) {
  i = pos[,1] == chr
  geno_chr = data.frame(ID = selected_individuals$V1,genotypes[,i])
  fwrite(geno_chr,file = file.path('genotypes_by_chromosome',sprintf('chr%02d.012.csv',chr)))
}

genotypes = as.matrix(genotypes)
pos$maf = colMeans(genotypes)/2
pos$maf[pos$maf>0.5] = 1-pos$maf[pos$maf>0.5]


genotypes_filtered = genotypes[,pos$maf > 0.05]
pos_filtered = pos[pos$maf > 0.05,]

LOO_ks = list()
try(dir.create('LOO_Ks'))


Ks = list()
for(chr in 1:10) {
  print(chr)
  mat = genotypes_filtered[,pos_filtered[,1] == chr]
  p = colMeans(mat,na.rm=T)/2
  mat = sweep(mat,2,2*p,'-')
  Ks[[chr]] = list(K = tcrossprod(mat)/ncol(mat),m = ncol(mat),c = 2*sum(p*(1-p)))
}


for(chr in 1:10) {
  print(chr)
  K = 0
  m = 0
  c = 0
  for(chr2 in 1:10) {
    if(chr == chr2) next
    Ki = Ks[[chr2]]
    K = K + Ki$K*Ki$m
    m = m + Ki$m
    c = c + Ki$c
  }
  # K = K/m
  K = K/c
  rownames(K) = colnames(K) = selected_individuals$V1
  LOO_ks[[chr]] = K
  write.csv(K,file = sprintf('LOO_Ks/K_chr_%02d.csv',chr))
}

K = 0
m = 0
c = 0
for(chr2 in 1:10) {
  Ki = Ks[[chr2]]
  K = K + Ki$K*Ki$m
  m = m + Ki$m
  c = c + Ki$c
}
# K = K/m
K = K/c
rownames(K) = colnames(K) = selected_individuals$V1
write.csv(K,file = sprintf('K_allChr.csv',chr))
