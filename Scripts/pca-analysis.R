####################################################################################
# PCA Analysis
#
# Author: Forrest Li
# Describe genetic variation in this population looking at PCs
####################################################################################

source(here::here('config.R'))
library(ggfortify)
# 
# K = fread('Genetic_data/Imputed_V4/K_allChr.csv', data.table = F)
# rownames(K) = K[,1]
# K = as.matrix(K[,-1])
cimmyt_grow = read.csv(here(env_data_dir, 'GEA-climate-nontransformed.csv'))
# 
# K = K[cimmyt_grow$Unique.ID, cimmyt_grow$Unique.ID]
# Kpca = prcomp(K)


genos = fread( 'Genetic_data/sampled_SNPs_forPCA_noNA.csv', header = T, data.table = F)
rownames(genos) = genos[,1]
sampled_genos = genos[,sample(ncol(genos), size = 5000 , replace = F)]
geno_pca = prcomp(sampled_genos, rank. = 10)
autoplot( geno_pca, data = cimmyt_grow, colour = 'elevation')



(geno_biplot = autoplot( geno_pca, data = cimmyt_grow, color = 'tmax'))
png(here(plot_dir, 'Manuscript', 'geno_biplot.png'), width = 491, height = 491)
print(geno_biplot)
dev.off()

var_explained = geno_pca$sdev^2 / sum(geno_pca$sdev^2)
(scree_plot = qplot(c(1:10), var_explained[1:10]) +
  geom_line() +
  xlab('Principal component') +
  ylab('Genetic variance explained') +
  ggtitle('Scree plot of population accessions') +
  ylim(0, max(var_explained)))
png(here(plot_dir, 'Manuscript', 'scree_plot.png'), width = 491, height = 491)
print(scree_plot)
dev.off()

sum(var_explained[1:10])

library(gridGraphics)
library(grid)

cimmyt_pca = cbind(cimmyt_grow, geno_pca$x[,1:5])
colnames(cimmyt_pca)[149:153] = c('PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5')
pca_cor_matrix= cor(cimmyt_pca[,c('PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5')], cimmyt_pca[,c('LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation')])
if(!file.exists('Analyses/Tables/pca_cor_matrix.csv')) {
  write.csv(pca_cor_matrix, 'Analyses/Tables/pca_cor_matrix.csv')
}

pca1_vs_elevation = ggplot(cimmyt_pca, aes(x = PCA1, y = elevation)) +
  geom_point() +
  ggtitle('PCA1 vs elevation')
pca1_vs_latitude = ggplot(cimmyt_pca, aes(x = PCA1, y = LatNew)) +
  geom_point() + 
  ggtitle('PCA1 vs latitude')

(supp4 = plot_grid(pca1_vs_elevation, pca1_vs_latitude, labels = 'AUTO', nrow = 2))
# figure1 = plot_grid(elevation_map, environmental_correlations, labels = 'AUTO', col = 1, row = 2)
png(here(plot_dir, 'Manuscript', 'supp4.png'), height = 982, width = 491)
dev.off()


ggpairs(cimmyt_pca[,c('LatNew', 'tmin', 'tmax', 'trange', 'precipTot', 'aridityMean', 'rhMean', 'elevation', 'PCA1', 'PCA2', 'PCA3', 'PCA4','PCA5')],
        columns = 9:13)
