library(brms)
library(data.table)
outdir = 'Resuls'
try(dir.create(outdir,showWarnings = FALSE))

data = fread('../../Phenotype_data/blups_deregressed_CV.csv',data.table = F)
trial_info = read.csv('../../Phenotype_data/Trial_info.csv')

dir = '/group/runciegrp2/Projects/SeeD/'
env_data = fread(paste0(dir, 'Env_data/GEA-climate-nontransformed.csv'), data.table = F)
geno_info = read.csv(file.path(dir, 'Phenotype_data/selected_genotypeIDs.csv'))
passport_data = readxl::read_excel(file.path(dir, 'Env_data/Passport/Maize climate data nov 2022 AEZ elevation by growing season and flowering period.xlsx'))
dna_to_tc_gid = as.data.frame(readxl::read_xlsx(file.path(dir, 'Phenotype_data/Hearne_data/dataset/Blups_01.xlsx'),sheet = 'DNA to TC GID'))


geno_info$AccID = dna_to_tc_gid[match(geno_info$Sample,dna_to_tc_gid$`Sample ID`),2]
geno_info$LatNew = passport_data$LatNew[match(geno_info$AccID,passport_data$GID)]
geno_info$LongNew = passport_data$LongNew[match(geno_info$AccID,passport_data$GID)]
geno_info$AltM = passport_data$AltM[match(geno_info$AccID,passport_data$GID)]

# geno_info = geno_info[!is.na(geno_info$LatNew) & geno_info$V1 %in% colnames(K) & geno_info$V1 %in% env_data$Unique.ID,]
# K = K[geno_info$V1,geno_info$V1]
env_data = env_data[match(geno_info$V1,env_data$Unique.ID),]

geno_info = merge(geno_info,env_data,by.x = 'V1',by.y = 'Unique.ID')


geno_env_variable = 'elevation'
trial_env_variable = 'Trial_elevation'
trait = 'GrainWeightPerHectareCorrected'
geno_env_variable = c('elevation','precipTot')[as.numeric(commandArgs(t=T)[1])]
trial_env_variable = c('Trial_elevation','precipTot')[as.numeric(commandArgs(t=T)[2])]
trait = unique(data$Trait)[as.numeric(commandArgs(t=T)[3])]


data_trait = subset(data,Trait==trait)
data_trait$geno_x = geno_info[[geno_env_variable]][match(data_trait$SampleID,geno_info$Sample)]
data_trait$trial_x = trial_info[[trial_env_variable]][match(data_trait$Experimento,trial_info$Experimento)]
data_trait = subset(data_trait,!is.na(geno_x))

std_x_mean = mean(data_trait$geno_x)
std_x_sd = sd(data_trait$geno_x)

data_trait$geno_x_std = (data_trait$geno_x - std_x_mean)/std_x_sd
data_trait$trial_x_std = (data_trait$trial_x - std_x_mean)/std_x_sd
data_trait$y = scale(data_trait$Value)

data_trait$ExpTester = interaction(data_trait$Experimento,data_trait$Tester,drop=TRUE)

prior1 = prior(normal(0,1),nlpar="c") + prior(normal(0,1),nlpar="a",ub=0) + prior(normal(0,3),nlpar="h")
bm = brm(bf(y ~ c+a*(geno_x_std-h)^2, a~ 0+Experimento,h~1+trial_x_std+(1|Experimento),sigma~0+Experimento,c~1+(1|ExpTester), nl = TRUE),
         data = data_trait, prior = prior1,
         control = list(adapt_delta = 0.9),cores=4)

saveRDS(list(bm=bm,data_trait=data_trait,std_x_mean=std_x_mean,std_x_sd=std_x_sd),file = file.path(outdir,sprintf('bm_%s_%s.rds',geno_env_variable,trait)))

# library(brms)
# library(ggplot2)
# files = list.files(path = outdir,pattern='.rds')
# pdf('Optimum_posteriors.pdf')
# for(file in files) {
#   trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
#   env = strsplit(file,'_')[[1]][[2]]
#   bm_list = readRDS(file.path(outdir,file))
#   bm = bm_list$bm
#   d = bm_list$data_trait
#   new_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
#   
#   
#   h_posteriors = as_draws_df(posterior_linpred(bm,newdata = new_data,nlpar='h'))
#   colnames(h_posteriors)[1:nrow(new_data)] = new_data$Experimento
#   h_posteriors_tall = tidyr::pivot_longer(h_posteriors[,1:nrow(new_data)], cols = 1:nrow(new_data))
#   h_posteriors_tall$trial_x = new_data$trial_x_std[match(h_posteriors_tall$name,new_data$Experimento)]
#   # h_posteriors_tall$geno_x = bm_list$data_trait$geno_x_std[match(h_posteriors_tall$name,bm_list$data_trait$SampleID)]
#   
#   lines = as_draws_df(bm,c('b_h_Intercept','b_h_trial_x_std'))
#   X = seq(min(bm_list$data_trait$trial_x_std),max(bm_list$data_trait$trial_x_std),length=100)
#   fit_df = do.call(rbind,lapply(X,function(x) {
#     y = unlist(lines[,1] + lines[,2]*x)
#     data.frame(x=x,low = quantile(y,.025),high=quantile(y,0.975),mean = mean(y),value=NA)
#   }))
#   
#   rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
#   
#   h_posteriors_tall$value = rescale(h_posteriors_tall$value,bm_list)
#   h_posteriors_tall$trial_x = rescale(h_posteriors_tall$trial_x,bm_list)
#   fit_df$x = rescale(fit_df$x,bm_list)
#   fit_df$low = rescale(fit_df$low,bm_list)
#   fit_df$high = rescale(fit_df$high,bm_list)
#   fit_df$mean = rescale(fit_df$mean,bm_list)
#   h_posteriors_summary = aggregate(value~trial_x+name,h_posteriors_tall,FUN=function(y) {
#     # browser()
#     c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
#   })
#   
#   h_posteriors_summary = data.frame(h_posteriors_summary[,-ncol(h_posteriors_summary)],as.data.frame(h_posteriors_summary[,ncol(h_posteriors_summary)]))
#   (p = ggplot() + ggtitle(trait) + xlab(env) +
#       geom_ribbon(data = fit_df,aes(x=x,ymin = low,ymax = high),alpha = 0.3)+
#       geom_line(data = fit_df,aes(x=x,y=mean))+
#       geom_abline(slope=1,intercept=0,color='red')+
#       geom_segment(data=h_posteriors_summary,aes(x=trial_x,y=`low.2.5.`,xend=trial_x,yend=`high.97.5.`,group=name,color=name))
#     # geom_jitter(alpha = 0.01,aes(y=value,color=name)) +
#     # geom_violin(aes(y=value,group = name,color=name),position = position_dodge())
#   );print(p)
#   
# }
# dev.off()
# 
# 
# library(brms)
# library(ggplot2)
# files = list.files(path = outdir,pattern='.rds')
# pdf('Experimento_posteriors.pdf')
# for(file in files) {
#   trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
#   env = strsplit(file,'_')[[1]][[2]]
#   bm_list = readRDS(file.path(outdir,file))
#   bm = bm_list$bm
#   d = bm_list$data_trait
#   trial_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
#   trial_data = trial_data[order(trial_data$trial_x),]
#   ribbons = c()
#   points = c()
#   for(exp in trial_data$Experimento) {
#     d_exp = subset(d,Experimento==exp)
#     X = seq(min(d_exp$geno_x_std),max(d_exp$geno_x_std),length=100)
#     newdata=data.frame(
#       Experimento=exp, ExpTester = d_exp$ExpTester[1],
#       geno_x_std = X,
#       trial_x_std = d_exp$trial_x_std[1]
#     )
#     newdata$ExpTester[-1] = NA
#     pre = posterior_epred(bm,newdata = newdata)
#     pre_summary = apply(pre,2,function(y) {
#       c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
#     })
#     pre_summary = data.frame(geno_x = X,t(pre_summary))
#     rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
#     pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
#     # pre_summary$low.2.5. = rescale(pre_summary$low.2.5.,bm_list)
#     # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
#     # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
#     ribbons = rbind(ribbons,data.frame(Experimento = exp,pre_summary))
#     points = rbind(points,data.frame(Experimento = exp,y = d_exp$y,geno_x = d_exp$geno_x))
#   }
#   ribbons$Experimento = factor(ribbons$Experimento,levels = trial_data$Experimento)
#   p=ggplot(ribbons) + ggtitle(trait) + xlab(env) + 
#     geom_ribbon(aes(x=geno_x,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
#     geom_point(data = points,aes(x = geno_x,y=y)) + 
#     geom_line(aes(x=geno_x,y=mean),linewidth=2,color='blue') +
#     facet_wrap(~Experimento,scales = 'free_y')
#   print(p)
# }
# dev.off()
