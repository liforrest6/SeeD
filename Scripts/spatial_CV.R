library(data.table)
library(spatialsample)
library(sf)
library(SpatialKDE)
library(sp)
library(dplyr)
library(tmap)
library(mgcv)
library(gamm4)
library(sommer)
library(lme4qtl)

K = fread('Genetic_data/Imputed_V4/K_allChr.csv',data.table=F)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
env_data = read.csv(here(env_data_dir, 'GEA-climate-invnormtransformed.csv'))
pheno_data = read.csv(here(phenotype_data_dir, 'blups_std.csv'))
geno_info = read.csv(here(genetic_data_dir, 'selected_genotypeIDs.csv'))
passport_data = readxl::read_excel(here(env_data_dir, 'cimmyt_data.xlsx'))
dna_to_tc_gid = as.data.frame(readxl::read_xlsx(here(phenotype_data_dir, 'Blups_01.xlsx') ,sheet = 'DNA to TC GID'))

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


library(rrBLUP)
library(foreach)
library(doParallel)
library(ggplot2)
library(cowplot)
registerDoParallel(10)
traits = colnames(geno_info)[6:12]

geno_info[,traits] = scale(geno_info[,traits])

# spatial clustering
clusters_spatial = spatial_clustering_cv(geno_sf,v=10,buffer = 4e5)
# clusters_spatial = clusters
pdf('spatial_cluster_map.pdf')
autoplot(clusters)
dev.off()

res_spclust = foreach(i=1:nrow(clusters),.combine = rbind) %dopar% {
  training = analysis(clusters$splits[[i]])
  testing = assessment(clusters$splits[[i]])
  
  d = geno_info
  testing_index = match(testing$V1,d$V1)
  training_index = match(training$V1,d$V1)
  d[-training_index,traits] = NA
  foreach(trait=traits,.combine = rbind) %do% {
    print(c(i,trait))
    m = mixed.solve(d[[trait]],K=K)
    u = m$u + m$beta[1]
    data.frame(fold=i,trait=trait,r = cor(geno_info[[trait]][testing_index],u[testing_index]),rmse = sqrt(mean((geno_info[[trait]][testing_index]-u[testing_index])^2)))
  }
}
clusters = clusters_spatial
res_spclust_sommer = res_spclust
# res_spclust_sommer = foreach(i=1:nrow(clusters),.combine = rbind) %dopar% {
#   training = analysis(clusters$splits[[i]])
#   testing = assessment(clusters$splits[[i]])
#   
#   d = geno_info
#   testing_index = match(testing$V1,d$V1)
#   training_index = match(training$V1,d$V1)
#   d[-training_index,traits] = NA
#   # d = geno_info[training_index,]
#   foreach(trait=traits,.combine = rbind) %do% {
#     print(c(i,trait))
#     # Qt = t(chol(K[training_index,training_index]))
#     d$y = d[[trait]]
#     m = mmer(y~1,random = ~vsr(V1,Gu=K[d$V1,d$V1]) + spl2Da(LongNew,LatNew,nsegments=c(50,50),degree=c(3,3),
#                                                     penaltyord=c(3,3)),
#              # rcov=~spl2Da(LongNew,LatNew,nsegments=c(20,20),degree=c(3,3),
#              #                   penaltyord=c(3,3)),
#              data = d,dateWarning = F)
#     data.frame(fold=i,trait=trait,accuracy = cor(geno_info[[trait]][testing_index],m$U[[1]]$y[geno_info$V1[testing_index]]))
#   }
# }

geno_sf_acc = geno_sf
for(i in 1:nrow(clusters)) {
  testing = assessment(clusters$splits[[i]])
  testing_index = match(testing$V1,geno_sf_acc$V1)
  for(trait. in traits) {
    geno_sf_acc[[trait.]][testing_index] = subset(res_spclust,trait==trait. & fold == i)$accuracy
  }
}

pdf('spatial_cluster_CV.pdf')
ggplot(res_spclust,aes(x=trait,y=accuracy)) + geom_boxplot(aes(group=trait))
ggplot(res_spclust,aes(x=trait,y=accuracy)) + geom_line(aes(group=fold,color=factor(fold)))
# plot(geno_sf_acc,key.pos=1)
for(trait in traits) {
  p1=ggplot(data.frame(geno_info,value = geno_sf[[trait]])) + 
    geom_point(aes(x=LongNew,y=LatNew,color=value)) + ggtitle(trait)
  p2=ggplot(data.frame(geno_info,accuracy = geno_sf_acc[[trait]])) + 
    geom_point(aes(x=LongNew,y=LatNew,color=accuracy))# + ggtitle(trait)
  print(plot_grid(p1,p2,nrow=1))
}
dev.off()



# random clustering
clusters_random = rsample::vfold_cv(geno_sf,v=10)
clusters_random = spatial_buffer_vfold_cv(geno_sf, radius = NULL, buffer = NULL)
# clusters_random = clusters
# pdf('random_cluster_map.pdf')
# autoplot(clusters)
# dev.off()
res_random = foreach(i=1:nrow(clusters),.combine = rbind) %dopar% {
  training = analysis(clusters$splits[[i]])
  testing = assessment(clusters$splits[[i]])
  
  d = geno_info
  testing_index = match(testing$V1,d$V1)
  training_index = match(training$V1,d$V1)
  d[-training_index,traits] = NA
  d = d[training_index,]
  Ki = K[d$V1,d$V1]
  sKi = svd(Ki)
  D = Diagonal(nrow(d),sKi$d)
  rownames(D) = colnames(D) = d$V1
  d$int = t(sKi$u) %*% matrix(1,nrow(d))
  foreach(trait=traits,.combine = rbind) %do% {
    print(c(i,trait))
    d$y = d[[trait]]
    d$qty = t(sKi$u) %*% d$y
    m = relmatLmer(qty~0+int+(1|V1),d,relmat = list(V1=D))
    u0 = sKi$u %*% t(m@optinfo$relmat$relfac$V1) %*% ranef(m)$V1[,1]
    u1 = as.matrix(K[testing_index,training_index] %*% sKi$u %*% (1/sKi$d * t(sKi$u)) %*% u0)[,1] + fixef(m)[1]
    data.frame(fold=i,trait=trait,r = cor(geno_info[[trait]][testing_index],u1),rmse = sqrt(mean((geno_info[[trait]][testing_index]-u1)^2)))
    # m = mixed.solve(d[[trait]],K=K)
    # data.frame(fold=i,trait=trait,accuracy = cor(geno_info[[trait]][testing_index],m$u[testing_index]))
  }
}

res_combined = bind_rows(data.frame(Type = 'spatial_K',res_spclust,Fold = sprintf('Fold%02d',res_spclust$fold)),
                     data.frame(Type = 'spatial_K+S',res_spclust_sommer,Fold = sprintf('Fold%02d',res_spclust$fold)),
                     data.frame(Type = 'random_K',res_random,Fold = NA))

geno_info$Fold = NA
for(i in 1:nrow(clusters_spatial)) {
  testing = assessment(clusters_spatial$splits[[i]])
  testing_index = match(testing$V1,geno_info$V1)
  geno_info$Fold[testing_index] = i
}

res_spclust$sd_trait = NA
for(i in 1:nrow(res_spclust)) {
  res_spclust$sd_trait[i] = sd(subset(geno_info,Fold==res_spclust$fold[i])[[trait]])
}


pdf('Spatial_vs_random_CV.pdf')
autoplot(clusters_spatial)
ggplot(subset(res_combined,Type != 'spatial_K+S'),aes(x=trait,y=r)) +geom_boxplot(aes(group = interaction(trait,Type),color=Type)) +  expand_limits(y=1)
ggplot(subset(res_combined,Type != 'spatial_K+S'),aes(x=trait,y=rmse)) +geom_boxplot(aes(group = interaction(trait,Type),color=Type)) +  expand_limits(y=1)
ggplot(geno_info,aes(x=LongNew,y=LatNew)) + geom_text(aes(label= Fold,color=factor(Fold)))
ggplot(subset(res_combined,Type == 'spatial_K'),aes(x=trait,y=r)) +
  geom_point(data = res_random) + 
  # geom_boxplot(aes(group = interaction(trait,Type)),position = position_dodge(width=0.9))+
  # geom_point(aes(color=Fold),position = position_jitterdodge(dodge.width=0.9,jitter.width=.01)) + 
  # geom_point(aes(group = interaction(trait,Type),color=Fold),
  #            position = position_jitterdodge(dodge.width=0.9,jitter.width=.1),
  #            size=3) +  #,position = position_dodge(width=0.9)
  geom_line(aes(group = interaction(Fold,Type),color=Fold)) +
  geom_text(aes(label=fold,group = interaction(trait,Type),color=Fold),
            size=3) +  ggtitle('spatial_K') +
  expand_limits(y=1)
ggplot(subset(res_combined,Type == 'spatial_K'),aes(x=trait,y=rmse)) +
  # geom_boxplot(aes(group = interaction(trait,Type)),position = position_dodge(width=0.9))+
  # geom_point(aes(color=Fold),position = position_jitterdodge(dodge.width=0.9,jitter.width=.01)) + 
  # geom_point(aes(group = interaction(trait,Type),color=Fold),
  #            position = position_jitterdodge(dodge.width=0.9,jitter.width=.1),
  #            size=3) +  #,position = position_dodge(width=0.9)
  geom_line(aes(group = interaction(Fold,Type),color=Fold)) +
  geom_text(aes(label=fold,group = interaction(trait,Type),color=Fold),
            size=3) +  ggtitle('spatial_K') +
  expand_limits(y=2)
ggplot(res_spclust,aes(x=sd_trait,y=r)) + geom_text(aes(label=fold,color=factor(fold))) + facet_wrap(~trait)
ggplot(res_spclust,aes(x=sd_trait,y=rmse)) + geom_text(aes(label=fold,color=factor(fold))) + facet_wrap(~trait)
dev.off()

pdf('Spatial_vs_random_CV_with_K_S.pdf')
autoplot(clusters_spatial)
ggplot(res_combined,aes(x=trait,y=r)) +geom_boxplot(aes(group = interaction(trait,Type),color=Type)) +  expand_limits(y=1)
ggplot(geno_info,aes(x=LongNew,y=LatNew)) + geom_text(aes(label= Fold,color=factor(Fold)))
ggplot(subset(res_combined,Type == 'spatial_K'),aes(x=trait,y=r)) +
  # geom_boxplot(aes(group = interaction(trait,Type)),position = position_dodge(width=0.9))+
  # geom_point(aes(color=Fold),position = position_jitterdodge(dodge.width=0.9,jitter.width=.01)) + 
  # geom_point(aes(group = interaction(trait,Type),color=Fold),
  #            position = position_jitterdodge(dodge.width=0.9,jitter.width=.1),
  #            size=3) +  #,position = position_dodge(width=0.9)
  geom_line(aes(group = interaction(Fold,Type),color=Fold)) +
  geom_text(aes(label=fold,group = interaction(trait,Type),color=Fold),
             size=3) +  ggtitle('spatial_K') +
  expand_limits(y=1)
ggplot(subset(res_combined,Type == 'spatial_K+S'),aes(x=trait,y=r)) +
  geom_line(aes(group = interaction(Fold,Type),color=Fold)) +
  geom_text(aes(label=fold,group = interaction(trait,Type),color=Fold),
            size=3) + ggtitle('spatial_K+S') +
  expand_limits(y=1)
res_comp = res_spclust
res_comp$K_S_r = res_spclust_sommer$r
ggplot(res_comp,aes(x=r,y=K_S_r)) + geom_point() + geom_abline(slope=1,intercept=0) + facet_wrap(~fold)
ggplot(res_comp,aes(x=r,y=K_S_r)) + geom_point() + geom_abline(slope=1,intercept=0) + facet_wrap(~trait)
dev.off()


# 
# 
# geno_sf_acc = geno_sf
# for(i in 1:nrow(clusters)) {
#   testing = assessment(clusters$splits[[i]])
#   testing_index = match(testing$V1,geno_sf_acc$V1)
#   for(trait. in traits) {
#     geno_sf_acc[[trait.]][testing_index] = subset(res_random,trait==trait. & fold == i)$accuracy
#   }
# }
# 
# pdf('random_cluster_CV.pdf')
# ggplot(res_random,aes(x=trait,y=accuracy)) + geom_boxplot(aes(group=trait))
# ggplot(res_random,aes(x=trait,y=accuracy)) + geom_line(aes(group=fold,color=factor(fold)))
# # plot(geno_sf_acc,key.pos=1)
# for(trait in traits) {
#   p1=ggplot(data.frame(geno_info,value = geno_sf[[trait]])) + 
#     geom_point(aes(x=LongNew,y=LatNew,color=value)) + ggtitle(trait)
#   p2=ggplot(data.frame(geno_info,accuracy = geno_sf_acc[[trait]])) + 
#     geom_point(aes(x=LongNew,y=LatNew,color=accuracy))# + ggtitle(trait)
#   print(plot_grid(p1,p2,nrow=1))
# }
# dev.off()


pdf('train_test_splits.pdf')
for(i in 1:nrow(clusters_spatial)) {
  training = analysis(clusters_spatial$splits[[i]])
  testing = assessment(clusters_spatial$splits[[i]])
  testing_index = match(testing$V1,geno_info$V1)
  training_index = match(training$V1,geno_info$V1)
  
  d = rbind(data.frame(geno_info[training_index,],Type='Training'),
            data.frame(geno_info[testing_index,],Type='Test'))
  
  print(ggplot(d,aes(x=LongNew,y=LatNew)) + geom_point(aes(color=Type)) + ggtitle(sprintf('Fold%02d',i)))
}
dev.off()


{
  i = 1
  training = analysis(clusters_random$splits[[i]])
  testing = assessment(clusters_random$splits[[i]])
  testing_index = match(testing$V1,geno_info$V1)
  training_index = match(training$V1,geno_info$V1)
  
  d = rbind(data.frame(geno_info[training_index,],Type='Training'),
            data.frame(geno_info[testing_index,],Type='Test'))
  
  print(ggplot(d,aes(x=LongNew,y=LatNew)) + geom_point(aes(color=Type)) + ggtitle(sprintf('Fold%02d',i)))
}
