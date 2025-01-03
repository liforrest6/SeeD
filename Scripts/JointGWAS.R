library(vcfR)
library(data.table)
library(tidyr)
library(dplyr)
library(JointGWAS)

library(Matrix)
library(foreach)
library(doParallel)
registerDoParallel(RcppParallel::defaultNumThreads()-1)

genotype_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome'
LOO_K_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/LOO_Ks'

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

if(file.exists(sprintf('%s/JointGWAS_output/GWAS_results_chr_%02d.csv',results_folder,chr))) q()

data_tall = fread(file.path(results_folder,'prepped_data.csv'),data.table=F)

# remove section when cholL_Sigma_inv already built
SampleIDs = unique(data_tall$SampleID)

try(dir.create(file.path(results_folder,'cholL_Sigma_inv')))

K_file = sprintf('%s/cholL_Sigma_inv/cholL_Sigma_inv_chr%02d.txt',results_folder,chr)
if(!file.exists(K_file)) {
  # K = fread('../../Genetic_data/allLines/K_allChr.csv',data.table=F)
  K = fread(sprintf('%s/K_chr_%02d.csv',LOO_K_dir,chr),data.table=F)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  K = K[SampleIDs,SampleIDs]

  # make K for FOAM
  diag(K) = 1
  K = K/4

  # Build covariance matrix
  Experimentos = unique(data_tall$Experimento)
  G_hat = read.csv(sprintf('%s/G_hat.csv',results_folder))
  rownames(G_hat) = G_hat[,1]
  G_hat = as.matrix(G_hat[,-1])[Experimentos,Experimentos]
  R_hat = read.csv(sprintf('%s/R_hat.csv',results_folder))
  rownames(R_hat) = R_hat[,1]
  R_hat = as.matrix(R_hat[,-1])[Experimentos,Experimentos]
  P_hat = read.csv(sprintf('%s/P_hat.csv',results_folder))
  rownames(P_hat) = P_hat[,1]
  P_hat = as.matrix(P_hat[,-1])[Experimentos,Experimentos]

  cholL_Sigma_inv = make_cholL_Sigma_inv(data_tall,'y','SampleID','Experimento',list(list(Row = K,Column = G_hat),list(Column = R_hat)),sparse=TRUE)
  fwrite(list(cholL_Sigma_inv[lower.tri(cholL_Sigma_inv,diag=T)]),file = K_file)
} else {
  cholL_Sigma_inv = matrix(0,nrow(data_tall),nrow(data_tall))
  cholL_Sigma_inv[lower.tri(cholL_Sigma_inv,diag=T)] = fread(K_file)[[1]]
  cholL_Sigma_inv = as(cholL_Sigma_inv,'dgCMatrix')
}

mat = fread(sprintf('%s/chr%02d.012.csv',genotype_dir,chr),data.table=F)
rownames(mat) = mat[,1]
mat = as.matrix(mat[,-1])
mat = mat[SampleIDs,]

maf = .5 - abs(.5-colMeans(mat)/2)

select_cols = maf > 0.01

mat = mat[,select_cols]

map = fread(file.path(genotype_dir,'../V2_V4_mapping.csv'),data.table=F)
map = map[match(colnames(mat),map$V4),]
map$MAF = maf[select_cols]

results = EMMAX_ANOVA(formula=y~Experimento:Tester + X+X:Experimento,data_tall,mat,
                      'SampleID',cholL_Sigma_inv,RcppParallel::defaultNumThreads()-1,
                      MAC_filter = 10,MAF_filter=0.01,mean_Imputed_filter = NULL)
anova = cbind(map,results$anova)
anova$Joint_PValue = pbeta(apply(anova[,c('X..Pvalue','X.Experimento..Pvalue')],1,min),1,2)
beta_hats = results$beta_hats
SEs = results$SEs
rownames(anova) = rownames(beta_hats) = rownames(SEs) = map$V4
try(dir.create(file.path(results_folder,'JointGWAS_output')))
write.csv(anova,file = sprintf('%s/JointGWAS_output/GWAS_results_chr_%02d.csv',results_folder,chr),row.names=T,quote=F)
write.csv(beta_hats,file = sprintf('%s/JointGWAS_output/GWAS_beta_hats_chr_%02d.csv',results_folder,chr),row.names=T,quote=F)
write.csv(SEs,file = sprintf('%s/JointGWAS_output/GWAS_SEs_chr_%02d.csv',results_folder,chr),row.names=T,quote=F)


# library(data.table)
# library(qqman)
# pdf('Manhattan_plots_Joint_G_GxE.pdf')
# pdf('qqPlots_Joint_G_GxE.pdf')
# traits = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")
# for(trait in traits) {
#   for(dataset in c(1,3)) {
#     for(rep in c(1,2)) {
#     for(transform in 2) {
#       name = paste(trait,dataset,transform,rep,sep='_')
#       print(name)
#       files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/JointGWAS_output',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'GWAS_results',full.names=T)
#       # break
#       if(length(files) <= 1) next
#       results = do.call(dplyr::bind_rows,lapply(files,function(x) fread(x,data.table=F)))
#       results = results[order(results$CHROM,results$POS),]
#
#       if('X..Pvalue' %in% colnames(results)) {
#         results$p = pbeta(apply(results[,c('X..Pvalue','X.Experimento..Pvalue')],1,min),1,2)
#       } else {
#         results$p = results$`Experimento.X..Pvalue`
#       }
#
#       p_value = 'p'
#       # results$p = results[[p_value]]
#       results$SNP = results$V4
#       results_sig = subset(results,results[[p_value]] < .02)
#       dev.set(2)
#       manhattan(results_sig,chr = 'CHROM',bp = 'POS',p = 'p',main = name)
#       dev.set(3)
#       # qq(results$p,main = name)
#       p=qqunif.plot(na.omit(results$p),
#                     main=paste(name,0.01),
#                     conf.alpha=.05
#       )
#       print(p)
#
#       # results = subset(results,MAF>0.05)
#       # results_sig = subset(results,results[[p_value]] < .02)
#       # dev.set(2)
#       # manhattan(results_sig,chr = 'CHROM',bp = 'POS',p = p_value,main = name)
#       # dev.set(3)
#       # # qq(results[[p_value]],main = name)
#       # p=qqunif.plot(na.omit(results[[p_value]]),
#       #               main=paste(name,0.05),
#       #               conf.alpha=.05
#       # )
#       # print(p)
#
#       # files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/GEMMA_univariate',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'p_wald_chr',full.names=T)
#
#     }}
#   }
# }
# graphics.off()
#
#
# # library(data.table)
# # library(qqman)
# # pdf('Manhattan_plots_Joint.pdf')
# # pdf('qqPlots_Joint.pdf')
# # traits = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")
# # for(trait in traits) {
# #   for(dataset in 1) {
# #     for(rep in 1) {
# #     for(transform in 2) {
# #       name = paste(trait,dataset,transform,sep='_')
# #       print(name)
# #       files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/JointGWAS_output',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'GWAS_results',full.names=T)
# #       # break
# #       if(length(files) <= 1) next
# #       results = do.call(dplyr::bind_rows,lapply(files,function(x) fread(x,data.table=F)))
# #       results = results[order(results$CHROM,results$POS),]
# #
# #       p_value = 'Experimento.X..Pvalue'
# #       results$p = results[[p_value]]
# #       results$SNP = results$V4
# #       results_sig = subset(results,results[[p_value]] < .02)
# #       dev.set(2)
# #       manhattan(results_sig,chr = 'CHROM',bp = 'POS',p = 'p',main = name)
# #       dev.set(3)
# #       # qq(results$p,main = name)
# #       p=qqunif.plot(na.omit(results$p),
# #                     main=paste(name,0.01),
# #                     conf.alpha=.05
# #       )
# #       print(p)
# #
# #       results = subset(results,MAF>0.05)
# #       results_sig = subset(results,results[[p_value]] < .02)
# #       dev.set(2)
# #       manhattan(results_sig,chr = 'CHROM',bp = 'POS',p = p_value,main = name)
# #       dev.set(3)
# #       # qq(results[[p_value]],main = name)
# #       p=qqunif.plot(na.omit(results[[p_value]]),
# #                     main=paste(name,0.05),
# #                     conf.alpha=.05
# #       )
# #       print(p)
# #
# #       files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/GEMMA_univariate',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'p_wald_chr',full.names=T)
# #
# #     }}
# #   }
# # }
# # graphics.off()
# #
# # genotype_dir = '/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome'
# # map = fread(file.path(genotype_dir,'../V2_V4_mapping.csv'),data.table=F)
# #
# #
# # library(data.table)
# # library(qqman)
# # pdf('Manhattan_plots_GEMMA.pdf')
# # pdf('qqPlots_GEMMA.pdf')
# # traits = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")
# # for(trait in traits) {
# #   for(dataset in 1) {
# #     for(rep in 1) {
# #       for(transform in 1:2) {
# #         name = paste(trait,dataset,transform,sep='_')
# #         print(name)
# #         files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/GEMMA_univariate',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'p_wald_chr',full.names=T)
# #         # break
# #         if(length(files) <= 1) next
# #         results = do.call(dplyr::bind_rows,lapply(files,function(x) fread(x,data.table=F)))
# #         map$p = apply(results[,-1],1,function(x) {
# #           pbeta(min(x,na.rm=T),1,sum(!is.na(x)))
# #         })[match(map$V4,results$V1)]
# #         map$SNP = map$V4
# #         map = subset(map,!is.na(p))
# #
# #         dev.set(2)
# #         manhattan(subset(map,p < 0.2),chr = 'CHROM',bp = 'POS',p = 'p',main = name)
# #         dev.set(3)
# #         qq(map$p,main = name)
# #
# #       }}
# #   }
# # }
# # graphics.off()
# #
# # trait = 'GrainWeightPerHectare'
# # # trait = 'DaysToFlowering'
# # # trait = 'BareCobWeight'
# # dataset=1
# # transform=2
# # rep=1
# # name = paste(trait,dataset,transform,sep='_')
# # print(name)
# # files = list.files(path = sprintf('JointGWAS/%s/dataset_%02d_%s/rep_%02d/GEMMA_univariate',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'p_wald_chr',full.names=T)
# # # break
# # if(length(files) <= 1) next
# # results = do.call(dplyr::bind_rows,lapply(files[!grepl('nullX',files)&!grepl('MVNperm',files)],function(x) fread(x,data.table=F)))
# # results_nullX = do.call(dplyr::bind_rows,lapply(files[grepl('nullX',files)],function(x) fread(x,data.table=F)))
# # results_MVNperm = do.call(dplyr::bind_rows,lapply(files[grepl('MVNperm',files)],function(x) fread(x,data.table=F)))
# #
# # traits = c('asi','dff','pcampo','pgrhah','polote','altpl','pgrha')
# # traitNames = c('ASI','DaysToFlowering','FieldWeight','GrainWeightPerHectareCorrected','BareCobWeight','PlantHeight','GrainWeightPerHectare')
# # H2s = as.data.frame(readxl::read_excel('../Phenotype_data/Burgueno/MASAGRO-Maize-Data.xlsx',sheet=3))
# # H2s$h2 = 1-H2s$`Average blup variance`/H2s$`Genetic variance`
# # H2s$Experimento = colnames(results)[match(H2s$Experiment,toupper(colnames(results)))]
# # H2s$TraitName = traitNames[match(H2s$Trait,traits)]
# # H2s = subset(H2s,!is.na(Experimento) & TraitName %in% traitNames)
# #
# #
# # qq_list = function(pvectors,truncation = 0.1,...) {
# #   # recover()
# #   qqs = lapply(pvectors,function(p) {
# #     p = na.omit(p)
# #     o = -log10(sort(p,decreasing=F))
# #     e = -log10(ppoints(length(p)))
# #     cbind(e=e,o=o)[e>-log10(truncation),]
# #   })
# #   xlim = c(0,max(do.call(c,lapply(qqs,function(x) x[,1]))))
# #   ylim = c(0,max(do.call(c,lapply(qqs,function(x) x[,2]))))
# #   def_args <- list(pch = 20, xlim = xlim, ylim = ylim, xlab = expression(Expected ~ ~-log[10](italic(p))),
# #                    ylab = expression(Observed ~ ~-log[10](italic(p))))
# #   dotargs <- list(...)
# #   tryCatch(do.call("plot", c(list(x = NA, y = NA), def_args[!names(def_args) %in%
# #                                                             names(dotargs)], dotargs)), warn = stop)
# #   sapply(1:length(qqs),function(i) points(qqs[[i]][,'e'],qqs[[i]][,'o'],col=i))
# #   legend('topleft',legend = names(pvectors),col = 1:length(qqs),pch=20)
# #   abline(0,1)
# # }
# #
# # h2s = H2s$h2[match(paste(colnames(results)[-1],trait),paste(H2s$Experimento,H2s$TraitName))]
# #
# # pdf(sprintf('GEMMA_%s_nullX_qq2.pdf',trait))
# # out=sapply(1+order(h2s),function(i) {
# #   # qq_list(list(X=results[,i],null_X = results_nullX[,i],MVNperm = results_MVNperm[,i]),main=paste(colnames(results)[i],H2s$h2[match(paste(colnames(results)[i],trait),paste(H2s$Experimento,H2s$TraitName))]))
# #   p=qqunif.plot(list(X=na.omit(results[,i]),null_X = na.omit(results_nullX[,i]),MVNperm = na.omit(results_MVNperm[,i])),
# #               main=paste(trait,colnames(results)[i],H2s$h2[match(paste(colnames(results)[i],trait),paste(H2s$Experimento,H2s$TraitName))]),
# #               conf.alpha=.05/length(h2s)
# #               )
# #   print(p)
# # })
# # dev.off()
# #
# # #https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
# library(lattice)
# qqunif.plot<-function(pvalues,
#                       should.thin=T, thin.obs.places=2, thin.exp.places=2,
#                       xlab=expression(paste("Expected (",-log[10], " p-value)")),
#                       ylab=expression(paste("Observed (",-log[10], " p-value)")),
#                       draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
#                       already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
#                       par.settings=list(superpose.symbol=list(pch=pch)), ...) {
#
#
#   #error checking
#   if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
#   if(!(class(pvalues)=="numeric" ||
#        (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
#     stop("pvalue vector is not numeric, can't draw plot")
#   if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
#   if (already.transformed==FALSE) {
#     if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
#   } else {
#     if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
#   }
#
#
#   grp<-NULL
#   n<-1
#   exp.x<-c()
#   if(is.list(pvalues)) {
#     nn<-sapply(pvalues, length)
#     rs<-cumsum(nn)
#     re<-rs-nn+1
#     n<-min(nn)
#     if (!is.null(names(pvalues))) {
#       grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
#       names(pvalues)<-NULL
#     } else {
#       grp=factor(rep(1:length(pvalues), nn))
#     }
#     pvo<-pvalues
#     pvalues<-numeric(sum(nn))
#     exp.x<-numeric(sum(nn))
#     for(i in 1:length(pvo)) {
#       if (!already.transformed) {
#         pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
#         exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
#       } else {
#         pvalues[rs[i]:re[i]] <- pvo[[i]]
#         exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
#       }
#     }
#   } else {
#     n <- length(pvalues)+1
#     if (!already.transformed) {
#       exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
#       pvalues <- -log10(pvalues)
#     } else {
#       exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
#     }
#   }
#
#
#   #this is a helper function to draw the confidence interval
#   panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
#     require(grid)
#     conf.points = min(conf.points, n-1);
#     mpts<-matrix(nrow=conf.points*2, ncol=2)
#     for(i in seq(from=1, to=conf.points)) {
#       mpts[i,1]<- -log10((i-.5)/n)
#       mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
#       mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
#       mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
#     }
#     grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
#   }
#
#   #reduce number of points to plot
#   if (should.thin==T) {
#     if (!is.null(grp)) {
#       thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
#                                 exp.x = round(exp.x, thin.exp.places),
#                                 grp=grp))
#       grp = thin$grp
#     } else {
#       thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
#                                 exp.x = round(exp.x, thin.exp.places)))
#     }
#     pvalues <- thin$pvalues
#     exp.x <- thin$exp.x
#   }
#   gc()
#
#   prepanel.qqunif= function(x,y,...) {
#     A = list()
#     A$xlim = range(x, y)*1.02
#     A$xlim[1]=0
#     A$ylim = A$xlim
#     return(A)
#   }
#
#   #draw the plot
#   xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
#          prepanel=prepanel, scales=list(axs="i"), pch=pch,
#          panel = function(x, y, ...) {
#            if (draw.conf) {
#              panel.qqconf(n, conf.points=conf.points,
#                           conf.col=conf.col, conf.alpha=conf.alpha)
#            };
#            panel.xyplot(x,y, ...);
#            panel.abline(0,1);
#          }, par.settings=par.settings, ...
#   )
# }
# #
# #
# # traits = c('asi','dff','pcampo','pgrhah','polote','altpl','pgrha')
# # traitNames = c('ASI','DaysToFlowering','FieldWeight','GrainWeightPerHectareCorrected','BareCobWeight','PlantHeight','GrainWeightPerHectare')
# # H2s = as.data.frame(readxl::read_excel('../Phenotype_data/Burgueno/MASAGRO-Maize-Data.xlsx',sheet=3))
# # H2s$h2 = 1-H2s$`Average blup variance`/H2s$`Genetic variance`
# # H2s$Experimento = colnames(results)[match(H2s$Experiment,toupper(colnames(results)))]
# # H2s$TraitName = traitNames[match(H2s$Trait,traits)]
# # H2s = subset(H2s,!is.na(Experimento) & TraitName %in% traitNames)
# #
# #
# # trait = 'GrainWeightPerHectare'
# # # trait = 'DaysToFlowering'
# # # trait = 'BareCobWeight'
# # pdf(sprintf('SKAT_unweighted_qq.pdf'))
# # for(trait in traitNames) {
# # dataset=1
# # transform=2
# # rep=1
# # name = paste(trait,dataset,transform,sep='_')
# # print(name)
# # files = list.files(path = sprintf('single_trait_SKAT/%s/dataset_%02d_%s/rep_%02d/',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'SKAT_unweighted_p_chr',full.names=T)
# # # break
# # if(length(files) <= 1) next
# # results = do.call(rbind,lapply(files[!grepl('nullX',files)&!grepl('MVNperm',files)],function(x) fread(x,data.table=F)))
# #
# # h2s = H2s$h2[match(paste(colnames(results)[-c(1:5)],trait),paste(H2s$Experimento,H2s$TraitName))]
# #
# # out=sapply(5+order(h2s),function(i) {
# #   # qq_list(list(X=results[,i],null_X = results_nullX[,i],MVNperm = results_MVNperm[,i]),main=paste(colnames(results)[i],H2s$h2[match(paste(colnames(results)[i],trait),paste(H2s$Experimento,H2s$TraitName))]))
# #   p=qqunif.plot(list(SKAT=na.omit(results[,i])),
# #                 main=paste(trait,colnames(results)[i],H2s$h2[match(paste(colnames(results)[i],trait),paste(H2s$Experimento,H2s$TraitName))]),
# #                 conf.alpha=.05/length(h2s)
# #   )
# #   print(p)
# # })
# # }
# # dev.off()
# #
# #
# # library(data.table)
# # library(qqman)
# # pdf('Manhattan_plots_unweighted_SKAT.pdf')
# # pdf('qqPlots_unweighted_SKAT.pdf')
# # traits = c("ASI","DaysToFlowering","FieldWeight","GrainWeightPerHectareCorrected","GrainWeightPerHectare","BareCobWeight","PlantHeight")
# # for(trait in traits) {
# #   for(dataset in 1) {
# #     for(rep in 1) {
# #       for(transform in 1:2) {
# #         name = paste(trait,dataset,transform,sep='_')
# #         print(name)
# #         files = list.files(path = sprintf('single_trait_SKAT/%s/dataset_%02d_%s/rep_%02d/',trait,dataset,ifelse(transform==1,'orig','INT'),rep),pattern = 'SKAT_unweighted_p_chr',full.names=T)
# #         # break
# #         if(length(files) <= 1) next
# #         results = do.call(rbind,lapply(files,function(x) fread(x,data.table=F)))
# #         results$p = apply(results[,-c(1:5)],1,function(x) {
# #           pbeta(min(x,na.rm=T),1,sum(!is.na(x)))
# #         })
# #         results = subset(results,!is.na(p) & p>0)
# #         results$POS = (results$start+results$end)/2
# #         results$SNP = results$Genes
# #
# #         dev.set(2)
# #         manhattan(subset(results,p < 0.2),chr = 'CHROM',bp = 'POS',p = 'p',main = name)
# #         dev.set(3)
# #         qq(results$p,main = name)
# #
# #       }}
# #   }
# # }
# # graphics.off()
#
#
# # # X_cov = model.matrix(~Experimento:Tester,data_tall)
# # # X_cov = X_cov[,colSums(X_cov^2)>0]
# # # X_base = model.matrix(~1+Experimento,data_tall)
# # # colnames(X_base) = sub('Experimento','',colnames(X_base))
# # # X_base = X_base[,colSums(X_base^2)>0]
# # # assign = c(1,rep(2,ncol(X_base)-1))
# #
# # # beta_hats = matrix(NA,nrow = ncol(mat),ncol = ncol(X_base),dimnames = list(map$snp,colnames(X_base)))
# # # SEs = matrix(NA,nrow = ncol(mat),ncol = ncol(X_base),dimnames = list(map$snp,colnames(X_base)))
# # # map$F = map$p = NA
# #
# # test_sets = rep(1:ncol(mat),each = ncol(X_base))
# # chunks = unique(c(seq(0,ncol(mat),by = 1000),ncol(mat)))
# # anova_results = c()
# # beta_hats = c()
# # SEs = c()
# # for(j in 2:length(chunks)) {
# #   index = seq(chunks[j-1]+1,chunks[j])
# #   print(range(index))
# #   results = fastEMMAX_ANOVAs(formula=y~Experimento:Tester + X+X:Experimento,data_tall,mat[,index],'SampleID',cholL_Sigma_inv,RcppParallel::defaultNumThreads()-1)
# #   anova_results = rbind(anova_results,data.frame(map[index,],results[[1]]$anova))
# #   beta_hats = rbind(beta_hats,results[[1]]$beta_hats)
# #   SEs = rbind(SEs,results[[1]]$SEs)
# #   X_test = c()
# #   X_test = do.call(cbind,lapply(index,function(i) X_base*mat[data_tall$SampleID,i]))
# #   test_sets = rep(1:length(index),each = ncol(X_base))
# #   results = fastEMMAX_F(data_tall$y,X_test,X_cov,chol_Sigma,test_sets,assign = assign,ML = F,REML=F,mc.cores = 1,verbose = F)
# #   gc()
# #   map$F[index] = results$Fvalues[,1]
# #   map$p[index] = results$Pvalues[,1]
# #   beta_hats[index,] = results$beta_hats
# #   SEs[index,] = results$SEs
# # }
# # write.csv(anova_results,file = sprintf('%s/GxE_GWAS_results_chr_%02d.csv',results_folder,chr))
# # write.csv(beta_hats,file = sprintf('%s/GxE_GWAS_beta_hats_chr_%02d.csv',results_folder,chr))
# # write.csv(SEs,file = sprintf('%s/GxE_GWAS_SEs_chr_%02d.csv',results_folder,chr))
#
#
# # library(data.table)
# # library(qqman)
# # pdf('Manhattan_plots_Joint.pdf')
# # pdf('qqPlots_Joint.pdf')
# # for(trait in colnames(data)[12:17]) {
# #   for(dataset in 1:3) {
# #     print(paste(trait,dataset,sep='_'))
# #     name = paste(trait,dataset,sep='_')
# #     files = list.files(path = 'Results',pattern = name,full.names=T)
# #     files = files[grep('Joint_GWAS',files)]
# #     if(length(files) <= 1) next
# #     results = do.call(bind_rows,lapply(files,function(x) fread(x,data.table=F)))
# #     results = subset(results,maf>0.05)
# #
# #     p_value = 'P'
# #     results_sig = subset(results,results[[p_value]] < .02)
# #     dev.set(2)
# #     manhattan(results_sig,chr = 'Chr',bp = 'pos',p = p_value,main = name)
# #     dev.set(3)
# #     qq(results[[p_value]],main = name)
# #   }
# # }
# # graphics.off()
# # #
# # #
# # # #
# # top_results = c()
# # for(trait in colnames(data)[12:17]) {
# #   for(dataset in 3) {
# #     print(paste(trait,dataset,sep='_'))
# #     name = paste(trait,dataset,sep='_')
# #     files = list.files(path = 'Results',pattern = name,full.names=T)
# #     files = files[grep('Joint_GWAS',files)]
# #     if(length(files) <= 1) next
# #     results = do.call(bind_rows,lapply(files,function(x) fread(x,data.table=F)))
# #     clumped_results = fread(sprintf('Results_clumped/Joint_GWAS_LOCO_results_%s_3_clumped.csv',trait),data.table=F)
# #     clumped_results = subset(clumped_results,P < 1e-5)
# #     results = results[match(clumped_results$SNP,results$snp),]
# #     results$Total_Clumped = clumped_results$TOTAL
# #     top_results = bind_rows(top_results,data.frame(Trait = trait,dataset=dataset,results))
# #   }
# # }
# # trials = read.csv('../Hearne_data/Trial_info.csv',stringsAsFactors = F)
# # trials = trials[order(trials$Cycle,trials$Trial_elevation),]
# # trials$Cycle = factor(trials$Cycle)
# # effect_plot = function(top_results,i) {
# #   trials_i = subset(trials,Experimento %in% colnames(top_results))
# #   trials_i$B = unlist(top_results[i,match(trials_i$Experimento,colnames(top_results))])
# #   trials_i$SE = unlist(top_results[i,match(paste0('SE_',trials_i$Experimento),colnames(top_results))])
# #   par(mar=c(7,2,2,1)+.1)
# #   bp=barplot(trials_i$B,beside=T,ylim = c(-1,1)*max(abs(trials_i$B)+2*abs(trials_i$SE),na.rm=T),
# #           main = paste(top_results$Trait[i],top_results$dataset[i],top_results$snp[i]),
# #           names.arg = trials_i$Experimento,
# #           col = trials_i$Cycle,
# #           las=2)
# #   abline(h=0)
# #   arrows(bp[,1],trials_i$B-2*trials_i$SE,bp[,1],trials_i$B+2*trials_i$SE,angle=90,length=.1,code=3)
# # }
# #
# # pdf('trial_effects.pdf')
# # for(trait in unique(top_results$Trait)) {
# #   top_trait = subset(top_results,Trait == trait & maf > 0.05)
# #   top_trait = top_trait[order(top_trait$P),]
# #   for(i in 1:min(10,nrow(top_trait))) {
# #     effect_plot(top_trait,i)
# #   }
# # }
# # # trait = 'DaysToAnthesis'
# # # top_trait = subset(top_results,Trait == trait & Chr == 4 & maf > 0.05)
# # # top_trait = top_trait[order(top_trait$P),]
# # # for(i in 1:min(5,nrow(top_trait))) {
# # #   effect_plot(top_trait,i)
# # # }
# # dev.off()
#
# #
# # top_effects = as.matrix(top_results[,-c(1:13)]) + top_results$X
# # top_results = top_results[,1:13]
# # colnames(top_effects) = sub('EX_X.','',colnames(top_effects),fixed=T)
# #
# # top_results_MAC = c()
# # top_results_MAF = c()
# # for(chr in 1:10) {
# #   print(chr)
# #   mat = prep_genos(chr)
# #   for(trait in colnames(data)[12:17]) {
# #     # for(dataset. in 1:3) {
# #     index = which(with(top_results,Trait == trait & Chr == chr))
# #     for(i in index) {
# #       # aggregate(data[[trait]]~Experiment,data,FUN = length)
# #       data_trait = data[!is.na(data[[trait]]),]
# #       data_trait$X = 2-2*mat[data_trait$SampleID,top_results$snp[i]]
# #       mac = aggregate(X~Experimento,data_trait,FUN=function(x) sum(x))
# #       rownames(mac) = mac[,1]
# #       mac = t(mac[,-1,drop = F])
# #       top_results_MAC = bind_rows(top_results_MAC,data.frame(
# #         Snp = top_results$snp[i],mac))
# #       maf = aggregate(X~Experimento,data_trait,FUN=function(x) mean(x)/2)
# #       rownames(maf) = maf[,1]
# #       maf = t(maf[,-1,drop = F])
# #       top_results_MAF = bind_rows(top_results_MAF,data.frame(
# #         Snp = top_results$snp[i],maf))
# #     }
# #   }
# # }
# # top_results_MAC = top_results_MAC[match(top_results$snp,top_results_MAC$Snp),]
# # top_results_MAF = top_results_MAF[match(top_results$snp,top_results_MAF$Snp),]
# #
# # trials = read.csv('../Hearne_data/Trial_info.csv',stringsAsFactors = F)
# # trials = trials[match(colnames(top_effects),trials$Experimento),]
# # x = which(with(top_results,P_GxE < 1e-5 & maf > 0.05))
# # pdf('top_effects.pdf')
# # for(i in x) {
# #   trials$b = top_effects[i,]
# #   trials$MAC = unlist(top_results_MAC[i,trials$Experimento])
# #   p = ggplot(trials,aes(x=Trial_elevation,y=b)) + geom_point(aes(color = Cycle,size = MAC)) +
# #     ggtitle(paste(top_results$Trait[i],top_results$dataset[i],top_results$snp[i])) +
# #     geom_hline(yintercept=0)
# #   print(p)
# # }
# # dev.off()
# #
#
