# script to make transfer distance plots

library(ggplot2)
library(cowplot)
library(data.table)
library(MASS)


# data = fread('../prepped_data/blups_std.csv',data.table=F)
# data = fread('../prepped_data/blups_deregressed.csv',data.table=F)
data = fread(file.path(phenotype_data_dir, 'blups_deregressed_CV.csv'),data.table=F)
# data = fread(file.path(phenotype_data_dir, 'blups_std.csv'),data.table=F)
# data = fread(file.path(phenotype_data_dir, 'blups_deregressed.csv'),data.table=F)
TC_info = fread(file.path(phenotype_data_dir, 'TestCross_passport.csv'),data.table=F)
sample_info = fread(file.path(phenotype_data_dir, 'SEEDGWAS_passport_environ_TC_info.csv'))
sample_info = sample_info[match(TC_info$SampleID,sample_info$SampleID),]
TC_info$altitude = sample_info$altitude*10

# TC_info = data.frame(latitude = sample_info$latitude.x,longitude = sample_info$longitude.x)
lats = na.omit(TC_info$latitude)
lats_fraction = fractions(lats-floor(lats))
lats_fraction = data.frame(row = which(!is.na(TC_info$latitude)),do.call(rbind,lapply(as.character(lats_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
lats_denominators = table(lats_fraction$X2)

longs = na.omit(TC_info$latitude)
longs_fraction = fractions(longs-floor(longs))
longs_fraction = data.frame(row = which(!is.na(TC_info$longitude)),do.call(rbind,lapply(as.character(longs_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
longs_denominators = table(longs_fraction$X2)

TC_info$Lat_Long_denominator[lats_fraction$row] = lats_fraction$X2

# TC_info = subset(TC_info,Lat_Long_denominator > 20)
trial_info = fread(file.path(phenotype_data_dir, 'Trial_info.csv'),data.table=F)
traitNames = unique(data$Trait)

regions = list(Mexico = "MEXICO",
               CentralAmerica = c(
                 "GUATEMALA","PANAMA","EL SALVADOR","NICARAGUA","COSTA RICA","HONDURAS"),
               SouthAmerica = c(
                 "BRAZIL","BOLIVIA","ARGENTINA","CHILE","GUYANA","FRENCH GUIANA","PARAGUAY","URUGUAY","VENEZUELA","COLOMBIA","SURINAME","ECUADOR","PERU"),
               Caribean = c(
                 "PUERTO RICO","DOMINICAN REPUBLIC","SAINT VINCENT AND THE GRENADINES","GRENADA","ANTIGUA AND BARBUDA","BARBADOS","CUBA","HAITI","VIRGIN ISLANDS (U.S.)","TRINIDAD AND TOBAGO","GUADELOUPE","MARTINIQUE","JAMAICA","VIRGIN ISLANDS (BRITISH)")
          )
regions_table = do.call(rbind,lapply(names(regions),function(reg) data.frame(Region = reg,Country = regions[[reg]])))

# data = subset(data,SampleID %in% subset(TC_info,Country == 'MEXICO')$SampleID)
data$Country_genotype = TC_info$Country[match(data$SampleID,TC_info$SampleID)]
data$Region_genotype = regions_table$Region[match(data$Country_genotype,regions_table$Country)]
data$Elevation_genotype = TC_info$elevation_germinate[match(data$SampleID,TC_info$SampleID)]
data$Elevation_genotype = TC_info$altitude[match(data$SampleID,TC_info$SampleID)]
data$Latitude_genotype = TC_info$latitude[match(data$SampleID,TC_info$SampleID)]
data$Env_genotype = data$Elevation_genotype
data$Env_trial = trial_info$Trial_elevation[match(data$Experimento,trial_info$Experimento)]
data$Env_distance = data$Env_genotype - data$Env_trial

# data = subset(data,Country_genotype == 'MEXICO')
# data = subset(data,Region_genotype %in% c('Mexico','CentralAmerica'))

traitNames = unique(data$Trait)
plots = list()
models = list()
for(trait in traitNames) {
  data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance))
  data_trait$y = resid(lm(Value~Experimento:Tester,data_trait)) * 100
  plots[[trait]] = ggplot(data_trait,aes(x=Env_distance,y = y, color = Env_trial)) +
    geom_smooth(aes(group = Experimento,color = Env_trial),se=F,method = 'lm',formula = y~poly(x,2)) +
    scale_color_continuous(name = 'Trial elevation', low = 'blue', high = 'red')+
    ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
    theme_bw() +
    theme(plot.title = element_text(size = 14),
          legend.title = element_text(size = 12)) +
    xlab('Transfer distance') +
    ylab('Trait residual')
  # plots[[trait]] = ggplot(data,aes(x=Elevation_difference,y = y)) + geom_smooth(aes(group = Experimento,color = Trial_elevation),se=F,span=1) + ggtitle(trait)
  # plots[[trait]] = ggplot(data,aes(x=Elevation_difference,y = y)) + geom_point(size=.3) + geom_smooth(se=F,method = 'lm',formula = y~poly(x,2)) + ggtitle(trait)
  # data_trait = subset(data,!is.na(y))
  models[[trait]] = data.frame(
    Trait = trait,
    Env = tapply(data_trait$Env_trial,data_trait$Experimento,mean),
    do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) coef(lm(y~poly(Env_distance,2),data_trait[i,]))))
  )
  colnames(models[[trait]])[-c(1:2)] = c('Intercept','X','X2')
}
# (transfer_plots = cowplot::plot_grid(plotlist = plots[3:6],nrow = 2, ncol=2))
(transfer_plots = ggarrange(plotlist = plots[3:6], common.legend = T, nrow = 2, ncol = 2, legend = 'right'))
# cowplot::save_plot(p,file = 'Fig_1_blups.pdf',base_height = 12)

(all_transfer_plots = cowplot::plot_grid(plotlist = plots[1:6],nrow = 3, ncol=2))
png(here(plot_dir, 'Manuscript', 'all_transfer_plots.png'), width = 750, height = 750)
print(all_transfer_plots)
dev.off()

## how many of these models have beta coefficient that is negative for local adaptation?
models

transfer_meta_analysis = function(model) {
  trial_count = nrow(model)
  trial_negative = nrow(model[model$X2 < 0,])
  return(c(trial_count, trial_negative/trial_count))
}

transfer_meta_analyses = lapply(models, transfer_meta_analysis)

# png('~/Documents/Projects/gea-adaptation/plots/fig5-transfer-plots.png', width = 982, height = 328)
# print(fig5)
# dev.off()

# models = do.call(rbind,models)
# ggplot(models,aes(x=Env,y=X)) + geom_point() + facet_wrap(~Trait,scales = 'free')
# ggplot(models,aes(x=Env/1000,y=X2)) + geom_point() + facet_wrap(~Trait,scales = 'free') + xlab('Elevation (100m) of trial') + ylab('Coefficient of (sample-trial elevation)^2')
# 
# 
# plots = list()
# models2 = list()
# for(trait in traitNames) {
#   data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance))
#   # data_trait$y = resid(lm(Value~Experimento:Tester,data_trait))
#   data_trait$y = data_trait$Value
#   data_trait$y_resid = resid(lm(Value~Experimento:Tester,data_trait))
# 
#   plots[[trait]] =
#   ggplot(data_trait,aes(x=Env_genotype,y = y_resid)) + geom_smooth(aes(group = Experimento,color = Env_trial),se=F,method = 'lm',formula = y~poly(x,2)) + ggtitle(trait)
#   # plots[[trait]] =
#   #   ggplot(data_trait,aes(x=Env_distance,y = y)) + geom_smooth(aes(group = Experimento,color = Env_trial),se=F,span=.75) + ggtitle(trait)
#   # plots[[trait]] = ggplot(data_trait,aes(x=Env_distance,y = y)) + geom_point(size=.3) + geom_smooth(se=F,method = 'lm',formula = y~poly(x,2)) + ggtitle(trait)
#   envs = tapply(data_trait$Env_trial,data_trait$Experimento,mean)
#   coefs = do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) {
#     t1 = data_trait$Tester[i][1]
#     m = lm(y~poly(Env_genotype,2)+Tester,data_trait[i,])#+Latitude_genotype
#     # m2 = lm(y~poly(Env_genotype,2):Tester+Latitude_genotype+Tester,data_trait[i,])#
#     x = do.call(seq,as.list(c(range(data_trait[i,]$Env_genotype),length=1000)))
#     # if(coef(m)[3] > 0) {
#     #   c(coef(m)[1:3],x[which.min(predict(m,newdata = data.frame(Env_genotype=x,Tester=t1,Latitude_genotype=0)))],unlist(anova(m)[1,4:5]))
#     # } else {
#       c(coef(m)[1:3],x[which.max(predict(m,newdata = data.frame(Env_genotype=x,Tester=t1,Latitude_genotype=0)))],unlist(anova(m)[1,4:5]))
#     # }
#   }))
#   coefs = as.data.frame(coefs)
#   colnames(coefs)[1:4] = c('Intercept','X','X2','zero')
#   models2[[trait]] = data.frame(
#     Trait = trait,
#     Experimento = rownames(coefs),
#     Season = trial_info$Cycle[match(rownames(coefs),trial_info$Experimento)],
#     Env = envs[rownames(coefs)],
#     coefs
#   )
# 
#   # coefs$Experimento = factor(rownames(coefs),levels = trials$Experimento[order(trials$Trial_elevation)])
# 
#   # ggplot(data_trait,aes(x=Elevation_difference,y = y))+
#   #   geom_vline(data = coefs,aes(xintercept = zero)) +
#   #   # geom_point(size=.1)+
#   #   geom_smooth(aes(group = Experimento,color = Trial_elevation),se=F,method = 'lm',formula = y~poly(x,2)) +
#   #   ggtitle(trait) + facet_grid(Experimento~.,scales = 'free_y') +
#   #   geom_hline(yintercept = 0,size=.2) #
# 
# 
# }
# (p = cowplot::plot_grid(plotlist = plots,ncol=2))
# 
# models2 = do.call(rbind,models2)
# (p = ggplot(models2,aes(x=Env,y=X2)) + geom_point(aes(color = Season)) + facet_wrap(~Trait,scales = 'free'))
# # cowplot::save_plot(p,file = 'Curvature_coefficients_CV.pdf')
# (p = ggplot(models2,aes(x=Env,y=zero)) + geom_point(aes(color = Season,size = -log10(Pr..F.))) + facet_wrap(~Trait,scales = 'free')) + geom_abline(slope=1,intercept=0)
# (p = ggplot(models2,aes(x=Env,y=F.value)) + geom_point(aes(color = Season)) + facet_wrap(~Trait))# + geom_abline(slope=1,intercept=0)
# (p = ggplot(models2,aes(x=Env,y=-log10(Pr..F.))) + geom_point(aes(color = Season)) + facet_wrap(~Trait))# + geom_abline(slope=1,intercept=0)
# library(qqman)
# qq(subset(models2,Trait=='FieldWeight')$Pr..F.)
# qq(subset(models2,Trait=='BareCobWeight')$Pr..F.)
# with(subset(models2,Trait=='FieldWeight'),summary(lm(zero~Env,weights=F.value)))
# with(subset(models2,Trait=='BareCobWeight'),summary(lm(zero~Env,weights=F.value)))
# # cowplot::save_plot(p,file = 'Curvature_zero.pdf')
# (p = ggplot(models2,aes(x=Env,y=zero)) + geom_point(aes(color = X2 > 0)) + facet_wrap(~Trait,scales = 'free')) + geom_abline(slope=1,intercept=0)
# # cowplot::save_plot(p,file = 'Env of_zero.pdf')
# (p = ggplot(subset(models2,Trait == 'GrainWeightPerHectare'),aes(x=Env,y=zero)) +
#     geom_point() + geom_abline(slope=1,intercept=0) +
#     xlab('Env of Trial') + ylab('Env of maximum') + theme_bw())
# 
# 
# trait = 'FieldWeight'
# data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance))
# # data_trait$y = resid(lm(Value~Experimento:Tester,data_trait))
# data_trait$y = data_trait$Value
# data_trait$y_resid = resid(lm(Value~Experimento:Tester,data_trait))
# data_trait = data_trait[order(data_trait$Env_trial),]
# data_trait$Experimento = factor(data_trait$Experimento,levels = unique(data_trait$Experimento))
# trial_info_trait = subset(trial_info,Experimento %in% levels(data_trait$Experimento))
# trial_info_trait$Experimento = factor(trial_info_trait$Experimento,levels(data_trait$Experimento))
# 
# 
# ggplot(subset(data_trait,Experimento == 'm2012BTORN'),aes(x=(Env_genotype),y = y_resid))  + xlab('Elevation of origin') + ylab('Field Weight residual (g)') +
#   geom_smooth(aes(group = Experimento,color = Env_trial),se=T,method = 'lm',formula = y~poly(x,2)) + #ggtitle(trait) +
#   geom_vline(data =  subset(trial_info_trait,Experimento == 'm2012BTORN'),aes(xintercept=(Trial_elevation)))+ facet_wrap(~Experimento) +
#   theme_cowplot() + background_grid(major = 'xy') + theme(legend.position = 'none')
# 
# ggplot(data_trait,aes(x=(Env_genotype),y = y_resid))  + xlab('Elevation of origin') + ylab('Field Weight residual (g)') +
#   geom_smooth(aes(group = Experimento,color = Env_trial),se=T,method = 'lm',formula = y~poly(x,2)) + #ggtitle(trait) +
#   geom_vline(data =  trial_info_trait,aes(xintercept=(Trial_elevation)))+ facet_wrap(~Experimento) +
#   theme_cowplot() + background_grid(major = 'xy') + theme(legend.position = 'none')
# 
# ggplot(data_trait,aes(x=(Env_genotype),y = y_resid))  +
#   geom_smooth(aes(group = Experimento,color = Env_trial),se=F,method = 'lm',formula = y~poly(x,2)) +
#   theme_cowplot() + background_grid(major = 'xy')
# 
# nlog10_trans = function() {
#   scales::trans_new('-log10',function(x) log10(x),function(x) 10^{-x},scales::log_breaks(10),domain = c(1e-100, Inf))
# }
# ggplot(subset(models2,Trait == trait),aes(x=Env,y=zero)) + geom_point(aes(size = -log10(Pr..F.))) + geom_abline(slope=1,intercept=0) +
#   scale_size_area() + theme_cowplot() + background_grid(major = 'xy') +
#   ylab(sprintf('Elevation of predicted maximum of %s',trait)) + xlab('Elevation of trial')
# 
# library(tidyverse)
# 
# map <- get_stamenmap( bbox = c(left = -115, bottom = -40, right = -35, top = 32), zoom = 4, maptype = "toner-lite")
# 
# ggmap(map)+
#   geom_point(data = TC_info,aes(x = longitude, y = latitude,fill=elevation_germinate),color="grey90",pch=21,alpha=0.5,size=2)+
#   scale_fill_continuous(low="green",high="blue",name = 'Elevation',limits = c(0,3100))
# 
# range(trial_info$Trial_latitude)
# range(trial_info$Trial_longitude)
# map_trial <- get_stamenmap( bbox = c(left = -115, bottom = 15, right = -90, top = 32), zoom = 5, maptype = "toner-lite")
# ggmap(map_trial)+
#   geom_point(data = trial_info,aes(x = Trial_longitude, y = Trial_latitude,fill=Trial_elevation),color="grey10",pch=21,alpha=0.5,size=5)+
#   scale_fill_continuous(low="green",high="blue",name = 'Elevation',limits = c(0,3100))
# 
# 
# plots = list()
# for(trait in traitNames) {
#   data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance))
#   data_trait$y = resid(lm(Value~Experimento:Tester,data_trait))
#   plots[[trait]] = ggplot(data_trait,aes(x=Latitude_genotype,y=y)) + geom_smooth(aes(color = Env_trial,group = Experimento),se=F) + ggtitle(trait)
#   # plots[[trait]] = ggplot(data_trait,aes(x=Elevation_genotype,y=y)) + geom_smooth(aes(color = Env_trial,group = Experimento),se=F) + ggtitle(trait)
# }
# (p = cowplot::plot_grid(plotlist = plots,ncol=2))
# 
# 
# trait_coefs = c()
# for(trait in traitNames) {
#   data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance) & Region_genotype %in% c('Mexico','SouthAmerica'))
#   data_trait$Value = data_trait$Value/sd(data_trait$Value)
#   X = model.matrix(~0+poly(Elevation_genotype,2),data_trait)
#   colnames(X) = c('slope','curve')
#   data_trait = cbind(data_trait,X)
#   for(trial in unique(data_trait$Experimento)) {
#     data_trait_trial = subset(data_trait,Experimento == trial & Region_genotype %in% names(regions)[1:3])
#     m1 = lm(Value~Tester + Region_genotype + (slope+curve):Region_genotype,data_trait_trial)
#     anova(m1)
#     coefs = summary(m1)$coef
#     coefs = coefs[grep(':',rownames(coefs)),]
#     Region = sapply(sub('Region_genotype','',rownames(coefs)),function(x) strsplit(x,':')[[1]][1])
#     value = sapply(rownames(coefs),function(x) strsplit(x,':')[[1]][2])
#     # degree = as.numeric(substr(sapply(rownames(coefs),function(x) strsplit(x,')',fixed=T)[[1]][2]),1,1))
#     trait_coefs = rbind(trait_coefs,data.frame(Trait = trait,Experimento = trial,Region=Region,value=value,Elev = trial_info$Trial_elevation[match(trial,trial_info$Experimento)],
#                                                coefs))
#   }
# }
# ggplot(subset(trait_coefs,value=='slope'),aes(x=Elev,y=Estimate)) + geom_point(aes(color=Region)) + facet_wrap(~Trait,scales = 'free') + geom_smooth(aes(color=Region),se=F,method='lm')
# ggplot(subset(trait_coefs,value=='curve'),aes(x=Elev,y=Estimate)) + geom_point(aes(color=Region)) + facet_wrap(~Trait,scales = 'free')
# # boxplot(trait_coefs$Estimate~trait_coefs$Trait,las=2)
