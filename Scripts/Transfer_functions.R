# script to make transfer distance plots

library(cowplot)
library(data.table)
library(MASS)
library(brms)
library(viridis)

source(here::here('config.R'))

# data = fread('../prepped_data/blups_std.csv',data.table=F)
# data = fread(file.path(phenotype_data_dir, 'blups_deregressed.csv'),data.table=F)
data = fread(file.path(phenotype_data_dir, 'blups_deregressed_CV.csv'),data.table=F)
# data = fread(file.path(phenotype_data_dir, 'blups_std.csv'),data.table=F)
cimmyt_grow = read.csv('Env_data/GEA-climate-nontransformed.csv')
worldclim_raw = read.csv('Env_data/SEEDGWAS_worldclim.csv')
non_dup_acc = read.csv(here(phenotype_data_dir, 'selected_genotypeIDs.csv'))

data = data %>% filter(Trial_Classification != 'Stress')

# TC_info = fread('Phenotype_data/SEEDGWAS_passport_environ_TC_info.csv')
# TC_info %>% pull(AccGID) %>% unique() %>% length()

TC_info = fread(file.path(phenotype_data_dir, 'TestCross_passport.csv'),data.table=F)

# using old testcross information that may include groups as well as collections
# {
#   sample_info = fread(file.path(phenotype_data_dir, 'SEEDGWAS_passport_environ_TC_info.csv'))
#   sample_info = sample_info[match(TC_info$SampleID,sample_info$SampleID),]
#   TC_info$altitude = sample_info$altitude*10
#   TC_info$meanTemp = sample_info$meanTemp
#   TC_info$annualPrecipitation = sample_info$annualPrecipitation * 10
#   TC_info$precipWarmQ = sample_info$precipWarmQ
#   TC_info$tempSeasonality = sample_info$tempSeasonality
#   TC_info$meanWarmTempQ = sample_info$meanWarmTempQ
# }

## using new cimmyt data
{
  TC_info = TC_info[match(TC_info$SampleID, worldclim_raw$Sample.ID),]
  sample_info = worldclim_raw[match(TC_info$SampleID, worldclim_raw$Sample.ID),]
  TC_info$altitude = cimmyt_grow[match(TC_info$SampleID, cimmyt_grow$Sample.ID),]$elevation
  TC_info$meanTemp = sample_info$meanTemp
  TC_info$annualPrecipitation = sample_info$annualPrecipitation
  TC_info$precipWarmQ = sample_info$precipWarmQ
  TC_info$tempSeasonality = sample_info$tempSeasonality
  TC_info$meanTempWarmQ = sample_info$meanTempWarmQ
  # data = data %>% filter(SampleID %in% cimmyt_grow$Sample.ID)
}

## coordinate selection
# {
#   # TC_info = data.frame(latitude = TC_info$latitude,longitude = TC_info$longitude)
#   lats = na.omit(TC_info$latitude)
#   lats_fraction = fractions(lats-floor(lats))
#   lats_fraction = data.frame(row = which(!is.na(TC_info$latitude)),do.call(rbind,lapply(as.character(lats_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
#   lats_denominators = table(lats_fraction$X2)
#   
#   longs = na.omit(TC_info$latitude)
#   longs_fraction = fractions(longs-floor(longs))
#   longs_fraction = data.frame(row = which(!is.na(TC_info$longitude)),do.call(rbind,lapply(as.character(longs_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
#   longs_denominators = table(longs_fraction$X2)
#   
#   TC_info$Lat_Long_denominator[lats_fraction$row] = lats_fraction$X2
#   TC_info = subset(TC_info,Lat_Long_denominator > 20)
# }

# TC_info = data.frame(latitude = sample_info$latitude.x,longitude = sample_info$longitude.x)
# lats = na.omit(TC_info$latitude)
# lats_fraction = fractions(lats-floor(lats))
# lats_fraction = data.frame(row = which(!is.na(TC_info$latitude)),do.call(rbind,lapply(as.character(lats_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
# lats_denominators = table(lats_fraction$X2)
# 
# longs = na.omit(TC_info$latitude)
# longs_fraction = fractions(longs-floor(longs))
# longs_fraction = data.frame(row = which(!is.na(TC_info$longitude)),do.call(rbind,lapply(as.character(longs_fraction),function(x) as.numeric(strsplit(x,'/')[[1]]))))
# longs_denominators = table(longs_fraction$X2)
# 
# TC_info$Lat_Long_denominator[lats_fraction$row] = lats_fraction$X2


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
# data$Elevation_genotype = TC_info$elevation_germinate[match(data$SampleID,TC_info$SampleID)]
data$Elevation_genotype = TC_info$altitude[match(data$SampleID,TC_info$SampleID)]
data$meanTemp_genotype = TC_info$meanTemp[match(data$SampleID,TC_info$SampleID)]
data$annualPrecipitation_genotype = TC_info$annualPrecipitation[match(data$SampleID,TC_info$SampleID)]
data$Latitude_genotype = TC_info$latitude[match(data$SampleID,TC_info$SampleID)]
data$meanTempWarmQ_genotype = TC_info$meanTempWarmQ[match(data$SampleID,TC_info$SampleID)]
data$precipWarmQ_genotype = TC_info$precipWarmQ[match(data$SampleID,TC_info$SampleID)]

trial_worldclim = read.csv('Phenotype_data/Trial_worldclim.csv')
trial_info = merge(trial_info, trial_worldclim, by = 'Experimento')
data$Elevation_trial = trial_info$Trial_elevation[match(data$Experimento,trial_info$Experimento)]
data$meanTemp_trial = trial_info$meanTemp[match(data$Experimento,trial_info$Experimento)]
data$annualPrecipitation_trial = trial_info$annualPrecipitation[match(data$Experimento,trial_info$Experimento)]
data$Latitude_genotype = TC_info$latitude[match(data$SampleID,TC_info$SampleID)]
data$meanTempWarmQ_trial = trial_info$meanTempWarmQ[match(data$Experimento,trial_info$Experimento)]
data$precipWarmQ_trial = trial_info$precipWarmQ[match(data$Experimento,trial_info$Experimento)]


## calculate distance in environmental distance
data$Elevation_distance = data$Elevation_genotype - data$Elevation_trial
data$meanTemp_distance = data$meanTemp_genotype - data$meanTemp_trial
data$annualPrecipitation_distance = data$annualPrecipitation_genotype - data$annualPrecipitation_trial
data$meanTempWarmQ_distance = data$meanTempWarmQ_genotype - data$meanTempWarmQ_trial
data$precipWarmQ_distance = data$precipWarmQ_genotype - data$precipWarmQ_trial


# data = subset(data,Country_genotype == 'MEXICO')
# data = subset(data,Region_genotype %in% c('Mexico','CentralAmerica'))

####################################################################################
## generate frequentist plots and models
####################################################################################

traitNames = unique(data$Trait)
elevation_plots = list()
elevation_models = list()
temp_plots = list()
temp_models = list()
precipitation_plots = list()
precipitation_models = list()


createTransferModels = function(climate_var, plot_formula, absolute_distance) {
  env_plots = list()
  env_models = list()
  model_data = data
  if(climate_var == 'elevation') {
    model_data$Env_genotype = data$Elevation_genotype
    model_data$Env_trial = data$Elevation_trial
    model_data$Env_distance = data$Elevation_distance
  } else if(climate_var == 'temperature') {
    model_data$Env_genotype = data$meanTemp_genotype
    model_data$Env_trial = data$meanTemp_trial
    model_data$Env_distance = data$meanTemp_distance
  } else if(climate_var == 'precipitation') {
    model_data$Env_genotype = data$annualPrecipitation_genotype
    model_data$Env_trial = data$annualPrecipitation_trial
    model_data$Env_distance = data$annualPrecipitation_distance
  }
  if(absolute_distance){
    model_data$Env_distance = abs(model_data$Env_distance)
  }
  # model_data = model_data %>% filter(!Experimento %in% c('m2012AOBES', 'm2012BTORR'))
  
  for(trait in traitNames) {
    ## plot transfer distance plots
    data_trait = subset(model_data %>% filter(SampleID %in% cimmyt_grow$Sample.ID),Trait==trait & !is.na(Value) & !is.na(Env_distance))
    data_trait$y = resid(lm(Value~Experimento:Tester,data_trait)) * 100
    env_plots[[trait]] = ggplot(data_trait,aes(x=Env_distance,y = y, color = Env_trial)) +
      geom_smooth(aes(group = Experimento,color = Env_trial),se=F,method = 'lm',formula = plot_formula) +
      scale_color_continuous(name = sprintf('Trial %s', climate_var), low = 'blue', high = 'red')+
      ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            legend.title = element_text(size = 10)) +
      xlab('Transfer distance') +
      ylab('Trait residual')
    
    if(plot_formula == 'y ~ poly(x, 2)')
    {
      env_coefs = as.data.frame(do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) {
        t1 = data_trait$Tester[i][1]
        m = lm(y~poly(Env_genotype,2)+Tester,data_trait[i,])
        # m2 = lm(y~poly(Env_genotype,2):Tester+Latitude_genotype+Tester,data_trait[i,])#
        x = do.call(seq,as.list(c(range(data_trait[i,]$Env_genotype),length=1000)))
        if(coef(m)[3] > 0) {
          c(coef(m)[1:3],x[which.min(predict(m,newdata = data.frame(Env_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
        } else {
          c(coef(m)[1:3],x[which.max(predict(m,newdata = data.frame(Env_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
        }
      })))
      colnames(env_coefs)[1:4] = c('Intercept','X','X2','vertex')
    } else if(plot_formula == 'y ~ x') {
      env_coefs = as.data.frame(do.call(rbind, lapply(unique(data_trait$Experimento),function(i) {
        m = lm(y~Env_distance,subset(data_trait, Experimento == i))
        coef(m)[2]
      })))
      colnames(env_coefs) = c('linear_slope')
    }
    env_models[[trait]] = data.frame(
      Trait = trait,
      Env = tapply(data_trait$Env_trial,data_trait$Experimento,mean),
      env_coefs
      # do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) coef(lm(y~poly(Elevation_distance,2),data_trait[i,]))))
    )
  }
  env_models = do.call(rbind, env_models)
  return(list(env_plots, env_models))
}

quadratic = formula(y~poly(x, 2))
linear = formula(y~x)
elevation_models_quadratic = createTransferModels('elevation', quadratic, absolute_distance = F)
temperature_models_quadratic = createTransferModels('temperature', quadratic, absolute_distance = F)
precipitation_models_quadratic = createTransferModels('precipitation', quadratic, absolute_distance =F)

elevation_models_linear = createTransferModels('elevation', linear, absolute_distance = T)
temperature_models_linear = createTransferModels('temperature', linear, absolute_distance = T)
precipitation_models_linear = createTransferModels('precipitation', linear, absolute_distance = T)

(elevation_transfer_plots = ggarrange(plotlist = elevation_models_quadratic[[1]][1:6], common.legend = T, nrow = 3, ncol = 2, legend = 'right'))
(temp_transfer_plots = ggarrange(plotlist = temperature_models_quadratic[[1]][1:6], common.legend = T, nrow = 3, ncol = 2, legend = 'right'))
(precipitation_transfer_plots = ggarrange(plotlist = precipitation_models_quadratic[[1]][1:6], common.legend = T, nrow = 3, ncol = 2, legend = 'right'))

(quadratic_plots = ggarrange(ggarrange(plotlist = elevation_models_quadratic[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             ggarrange(plotlist = temperature_models_quadratic[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             ggarrange(plotlist = precipitation_models_quadratic[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             common.legend = F, nrow = 3, labels = c('A', 'B', 'C')))
(linear_plots = ggarrange(ggarrange(plotlist = elevation_models_linear[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             ggarrange(plotlist = temperature_models_linear[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             ggarrange(plotlist = precipitation_models_linear[[1]][3:5], common.legend = T, nrow = 1, ncol = 3, legend = 'right'),
                             common.legend = F, nrow = 3, labels = c('A', 'B', 'C')))

# png(here(plot_dir, 'Manuscript', 'elevation_transfer_plots.png'), width = 700, height = 700)
# print(elevation_transfer_plots)
# dev.off()
# png(here(plot_dir, 'Manuscript', 'temp_transfer_plots.png'), width = 700, height = 700)
# print(temp_transfer_plots)
# dev.off()
# png(here(plot_dir, 'Manuscript', 'precipitation_transfer_plots.png'), width = 700, height = 700)
# print(precipitation_transfer_plots)
# dev.off()
# 
# png(here(plot_dir, 'Manuscript', 'quadratic_transfer_plots.png'), width = 700, height = 700)
# print(quadratic_plots)
# dev.off()
# png(here(plot_dir, 'Manuscript', 'linear_transfer_plots.png'), width = 700, height = 700)
# print(linear_plots)
# dev.off()

# (transfer_plots = cowplot::plot_grid(plotlist = plots[3:6],nrow = 2, ncol=2))
# (transfer_plots = ggarrange(plotlist = list(elevation_models_quadratic[[1]][3][[1]], 
#                                             elevation_models_quadratic[[1]][4][[1]] + ggtitle(label = 'Grain Weight per Hectare'),
#                                             elevation_models_quadratic[[1]][5][[1]]), 
#                             common.legend = T, nrow = 1, ncol = 3, legend = 'right'))


# cowplot::save_plot(p,file = 'Fig_1_blups.pdf',base_height = 12)

# (all_transfer_plots = cowplot::plot_grid(plotlist = plots[1:6],nrow = 3, ncol=2))
# png(here(plot_dir, 'Manuscript', 'all_transfer_plots.png'), width = 750, height = 750)
# print(all_transfer_plots)
# dev.off()

## how many of these models have beta coefficient that is negative for local adaptation?

transfer_meta_analysis = function(model) {
  trial_count = nrow(model)
  trial_negative = nrow(model[model$X2 < 0,])
  return(paste(trial_negative, trial_count, sep = '/'))
}

####################################################################################
## frequentist optimality plots
####################################################################################

yield_traits = c('FieldWeight', 'GrainWeightPerHectareCorrected','BareCobWeight')
## how many of these trials have apex at 0?  (real optimality)
(elevation_optimal_plots = ggplot(elevation_models_quadratic[[2]] %>% filter(Trait %in% yield_traits),
                                  aes(x=Env,y=vertex)) + 
    geom_point(aes(color = X2 > 0, size = -log10(Pr..F.))) + 
    facet_wrap(~Trait,scales = 'free') + 
    geom_abline(slope=1,intercept=0) +
    geom_smooth(method = 'lm', se = T, data = elevation_models_quadratic[[2]] %>% filter(Trait %in% yield_traits) %>% filter(X2 < 0)) +
    theme_bw() +
    xlab('Trial elevation') +
    ylab('Polynomial max/min elevation') +
    scale_color_discrete(name = 'Polynomial coeff sign', labels = c('Negative', 'Positive'))) 

(temp_optimal_plots = ggplot(temperature_models_quadratic[[2]] %>% filter(Trait %in% yield_traits),
                             aes(x=Env,y=vertex)) + 
    geom_point(aes(color = X2 > 10, size = -log10(Pr..F.))) + 
    facet_wrap(~Trait,scales = 'free') + 
    geom_abline(slope=1,intercept=0) +
    geom_smooth(method = 'lm', se = T, data = temperature_models_quadratic[[2]] %>% filter(Trait %in% yield_traits) %>% filter(X2 < 0)) +
    theme_bw() +
    xlab('Trial mean temperature') +
    ylab('Polynomial max/min temperature') +
    scale_color_discrete(name = 'Polynomial coeff sign', labels = c('Negative', 'Positive')))

(precipitation_optimal_plots = ggplot(precipitation_models_quadratic[[2]] %>% filter(Trait %in% yield_traits),
                                      aes(x=Env,y=vertex)) + 
    geom_point(aes(color = X2 > 0, size = -log10(Pr..F.))) + 
    facet_wrap(~Trait,scales = 'free') + 
    geom_abline(slope=1,intercept=0) +
    geom_smooth(method = 'lm', se = T) +
    theme_bw() +
    xlab('Trial annual precipitation') +
    ylab('Polynomial max/min precipitation') +
    scale_color_discrete(name = 'Polynomial coeff sign', labels = c('Negative', 'Positive')))

# (transfer_optimality = ggarrange(elevation_optimal_plots, temp_optimal_plots, common.legend = T, labels = c('A', 'B'), nrow = 2))
# png(here(plot_dir, 'Manuscript', 'transfer_optimality.png'), width = 700, height = 700)
# print(transfer_optimality)
# dev.off()


i = 0
env_variable = c('elevation', 'meanTemp', 'annualPrecipitation')
for(model in list(elevation_models_quadratic[[2]], temperature_models_quadratic[[2]], precipitation_models_quadratic[[2]])) {
  i = i+1
  for(trait in yield_traits) {
    print(env_variable[i])
    print(trait)
    print(transfer_meta_analysis(model %>% filter(Trait == eval(trait))))
    model_anova = with(subset(model,Trait==trait & X2 < 0),summary(lm(vertex~Env, weights = F.value)))
    # print(model_anova$coefficients[2,3])
    print(model_anova$coefficients[2,4])
    # print(transfer_meta_analysis(model %>% filter(Trait == eval(trait))))
    # with(subset(elevation_models,Trait==trait & X2 < 0),summary(lm(vertex~Env, weights = F.value)))
  }
}


test = with(subset(elevation_models,Trait=='FieldWeight' & X2 < 0),summary(lm(vertex~Env)))

## anova tests to test if in each trait, trials follow 1:1 with optimal vertex (zero transfer distance)
with(subset(elevation_models,Trait=='FieldWeight' & X2 < 0),summary(lm(vertex~Env)))
with(subset(elevation_models,Trait=='BareCobWeight' & X2 < 0),summary(lm(vertex~Env)))
with(subset(elevation_models,Trait=='GrainWeightPerHectareCorrected' & X2 < 0),summary(lm(vertex~Env)))
with(subset(temp_models,Trait=='FieldWeight' & X2 < 0),summary(lm(vertex~Env)))
with(subset(temp_models,Trait=='BareCobWeight' & X2 < 0),summary(lm(vertex~Env)))
with(subset(temp_models,Trait=='GrainWeightPerHectareCorrected' & X2 < 0),summary(lm(vertex~Env)))
with(subset(precipitation_models,Trait=='FieldWeight' & X2 < 0),summary(lm(vertex~Env,weights=F.value)))
with(subset(precipitation_models,Trait=='GrainWeightPerHectareCorrected' & X2 < 0),summary(lm(vertex~Env,weights=F.value)))

# m2011BCL
subset(elevation_models, Trait == 'FieldWeight' & X2 < 0)
ggplot(subset(data %>% filter(SampleID %in% cimmyt_grow$Sample.ID),
              Trait=='FieldWeight' & !is.na(Value) & !is.na(Elevation_distance) & Experimento == 'm2012AAFTU'),
       aes(x=Elevation_distance,y = Value, color = Elevation_trial)) +
  geom_point() +
  geom_smooth(aes(group = Experimento,color = Elevation_trial),se=F,method = 'lm',formula = y~poly(x,2)) +
  # scale_color_continuous(name = 'Trial elevation', low = 'blue', high = 'red')+
  # ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  xlab('Transfer distance') +
  ylab('Trait residual')

t.test(elevation_models %>% filter(Trait == 'FieldWeight') %>% pull(vertex))
with(subset(elevation_models,Trait=='FieldWeight' & X2 < 0),summary(lm(vertex~Env)))


####################################################################################
## generate brms models for transfer plots to model 
####################################################################################

climate_var = 'precipitation'
model_data = data
if(climate_var == 'elevation') {
  model_data$Env_genotype = data$Elevation_genotype
  model_data$Env_trial = data$Elevation_trial
  model_data$Env_distance = data$Elevation_distance
} else if(climate_var == 'temperature') {
  model_data$Env_genotype = data$meanTemp_genotype
  model_data$Env_trial = data$meanTemp_trial
  model_data$Env_distance = data$meanTemp_distance
} else if(climate_var == 'precipitation') {
  model_data$Env_genotype = data$annualPrecipitation_genotype
  model_data$Env_trial = data$annualPrecipitation_trial
  model_data$Env_distance = data$annualPrecipitation_distance
}

trait = 'PlantHeight'
{
  data_trait = subset(model_data,Trait==trait)
  # data_trait$geno_x = geno_info[[geno_env_variable]][match(data_trait$SampleID,geno_info$Sample)]
  # data_trait$trial_x = trial_info[[trial_env_variable]][match(data_trait$Experimento,trial_info$Experimento)]
  data_trait = data_trait %>% filter(!is.na(Env_genotype))
  
  std_x_mean = mean(data_trait$Env_genotype)
  std_x_sd = sd(data_trait$Env_genotype)
  
  data_trait$Env_genotype_std = (data_trait$Env_genotype - std_x_mean)/std_x_sd
  data_trait$Env_trial_std = (data_trait$Env_trial - std_x_mean)/std_x_sd
  data_trait$y = scale(data_trait$Value)
  
  data_trait$ExpTester = interaction(data_trait$Experimento,data_trait$Tester,drop=TRUE)
  
  prior1 = prior(normal(0,1),nlpar="c") + prior(normal(0,1),nlpar="a",ub=0) + prior(normal(0,3),nlpar="h")
  bm = brm(bf(y ~ c+a*(Env_genotype_std-h)^2, a~ 0+Experimento,h~1+Env_trial_std+(1|Experimento),sigma~0+Experimento,c~1+(1|ExpTester), nl = TRUE),
           data = data_trait, prior = prior1,
           control = list(adapt_delta = 0.9),cores=4)
  
  saveRDS(list(bm=bm,data_trait=data_trait,std_x_mean=std_x_mean,std_x_sd=std_x_sd),file = file.path(here(analyses_dir, 'TransferPlots', sprintf('bm_%s_%s.rds',climate_var,trait))))
  sprintf('done %s %s', trait, climate_var)
}

files = list.files(here(analyses_dir, 'TransferPlots'), pattern='.rds',full.names = T)
slope_results = list()
for(file in files) {
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  new_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  
  h_posteriors = as_draws_df(posterior_linpred(bm,newdata = new_data,nlpar='h'))
  colnames(h_posteriors)[1:nrow(new_data)] = new_data$Experimento
  h_posteriors_tall = tidyr::pivot_longer(h_posteriors[,1:nrow(new_data)], cols = 1:nrow(new_data))
  h_posteriors_tall$Env_trial = new_data$Env_trial_std[match(h_posteriors_tall$name,new_data$Experimento)]
  # h_posteriors_tall$geno_x = bm_list$data_trait$geno_x_std[match(h_posteriors_tall$name,bm_list$data_trait$SampleID)]
  
  lines = as_draws_df(bm,c('b_h_Intercept','b_h_Env_trial_std'))
  # slope_posteriors = rbind(slope_posteriors,data.frame(Trait = trait,b_h_x_trial = lines$b_h_Env_trial_std))
  this_slope_posteriors = data.frame(Trait = trait,Env=env,b_h_x_trial = lines$b_h_Env_trial_std)
  slope_posteriors = rbind(slope_posteriors, this_slope_posteriors)
  X = seq(min(bm_list$data_trait$Env_trial_std),max(bm_list$data_trait$Env_trial_std),length=100)
  fit_df = do.call(rbind,lapply(X,function(x) {
    y = unlist(lines[,1] + lines[,2]*x)
    data.frame(x=x,low = quantile(y,.025),high=quantile(y,0.975),mean = mean(y),value=NA)
  }))
  
  rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
  
  h_posteriors_tall$value = rescale(h_posteriors_tall$value,bm_list)
  h_posteriors_tall$Env_trial = rescale(h_posteriors_tall$Env_trial,bm_list)
  fit_df$x = rescale(fit_df$x,bm_list)
  fit_df$low = rescale(fit_df$low,bm_list)
  fit_df$high = rescale(fit_df$high,bm_list)
  fit_df$mean = rescale(fit_df$mean,bm_list)
  h_posteriors_summary = aggregate(value~Env_trial+name,h_posteriors_tall,FUN=function(y) {
    # browser()
    c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
  })
  
  h_posteriors_summary = data.frame(h_posteriors_summary[,-ncol(h_posteriors_summary)],as.data.frame(h_posteriors_summary[,ncol(h_posteriors_summary)]))
  h_posteriors_summary$trial_x = h_posteriors_summary$Env_trial + rnorm(nrow(h_posteriors_summary),0,30)
  trait_slope_model = lm(mean ~ Env_trial, h_posteriors_summary)
  slope.emt = emtrends(trait_slope_model, ~1, var='Env_trial')
  slope_value = summary(slope.emt)$Env_trial.trend
  slope_results[[file]] = data.frame(trait = trait,
             env = env,
             # correlation = cor(h_posteriors_summary$Env_trial, h_posteriors_summary$mean) ^ 2,
             # slope_value = summary(slope.emt)$Env_trial.trend,
             # lower.CL = summary(slope.emt)$lower.CL,
             # upper.CL = summary(slope.emt)$upper.CL,
             slope_posteriors_mean = mean(this_slope_posteriors$b_h_x_trial),
             slope_posteriors_low = quantile(this_slope_posteriors$b_h_x_trial,.025),
             slope_posteriors_high = quantile(this_slope_posteriors$b_h_x_trial,0.975),
             slope_posteriors_p_val = t.test(this_slope_posteriors$b_h_x_trial, mu = 1)$p.value
             )
}
slope_results_df = bind_rows(slope_results)
t.test(slope_posteriors$b_h_x_trial, mu = 1)$p.value

trait_slope_model = lm(mean ~ Env_trial, h_posteriors_summary)
slope.emt = emtrends(trait_slope_model, ~1, var='Env_trial')
slope.means = emmeans(trait_slope_model, specs = 'Env_trial')
slope.emt

####################################################################################
## plotting
####################################################################################
files = list.files(here(analyses_dir, 'TransferPlots'), pattern='.rds',full.names = T)
optimum_posteriors = list()
slope_posteriors = c()
for(file in files) {
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  new_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  
  h_posteriors = as_draws_df(posterior_linpred(bm,newdata = new_data,nlpar='h'))
  colnames(h_posteriors)[1:nrow(new_data)] = new_data$Experimento
  h_posteriors_tall = tidyr::pivot_longer(h_posteriors[,1:nrow(new_data)], cols = 1:nrow(new_data))
  h_posteriors_tall$Env_trial = new_data$Env_trial_std[match(h_posteriors_tall$name,new_data$Experimento)]
  # h_posteriors_tall$geno_x = bm_list$data_trait$geno_x_std[match(h_posteriors_tall$name,bm_list$data_trait$SampleID)]
  
  lines = as_draws_df(bm,c('b_h_Intercept','b_h_Env_trial_std'))
  # slope_posteriors = rbind(slope_posteriors,data.frame(Trait = trait,b_h_x_trial = lines$b_h_Env_trial_std))
  slope_posteriors = rbind(slope_posteriors,data.frame(Trait = trait,Env=env,b_h_x_trial = lines$b_h_Env_trial_std))
  X = seq(min(bm_list$data_trait$Env_trial_std),max(bm_list$data_trait$Env_trial_std),length=100)
  fit_df = do.call(rbind,lapply(X,function(x) {
    y = unlist(lines[,1] + lines[,2]*x)
    data.frame(x=x,low = quantile(y,.025),high=quantile(y,0.975),mean = mean(y),value=NA)
  }))
  
  rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
  
  h_posteriors_tall$value = rescale(h_posteriors_tall$value,bm_list)
  h_posteriors_tall$Env_trial = rescale(h_posteriors_tall$Env_trial,bm_list)
  fit_df$x = rescale(fit_df$x,bm_list)
  fit_df$low = rescale(fit_df$low,bm_list)
  fit_df$high = rescale(fit_df$high,bm_list)
  fit_df$mean = rescale(fit_df$mean,bm_list)
  h_posteriors_summary = aggregate(value~Env_trial+name,h_posteriors_tall,FUN=function(y) {
    # browser()
    c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
  })
  
  h_posteriors_summary = data.frame(h_posteriors_summary[,-ncol(h_posteriors_summary)],as.data.frame(h_posteriors_summary[,ncol(h_posteriors_summary)]))
  h_posteriors_summary$trial_x = h_posteriors_summary$Env_trial + rnorm(nrow(h_posteriors_summary),0,30)
  (p = ggplot() + ggtitle(sprintf('%s, %s', trait, env)) + 
      xlab(sprintf('Trial %s', env)) +
      ylab(sprintf('Optimum %s', env)) +
      geom_ribbon(data = fit_df,aes(x=x,ymin = low,ymax = high),alpha = 0.3)+
      geom_line(data = fit_df,aes(x=x,y=mean))+
      geom_abline(slope=1,intercept=0,color='red')+
      geom_segment(data=h_posteriors_summary,aes(x=Env_trial,y=`low.2.5.`,xend=Env_trial,yend=`high.97.5.`,group=name,color=name))+
      geom_point(data=h_posteriors_summary,aes(x=Env_trial,y=mean,color=name),size=2) +
      theme(text = element_text(size = 5)) +
      theme_bw() 
    # geom_jitter(alpha = 0.01,aes(y=value,color=name)) +
    # geom_violin(aes(y=value,group = name,color=name),position = position_dodge())
  )
  # print(p)
  optimum_posteriors[[file]] = p
}

(slope_posterior_plot = ggplot(slope_posteriors,aes(x=Trait,y=b_h_x_trial, color = Trait)) + 
  geom_boxplot(aes(group = Trait)) + 
  geom_hline(yintercept=0) + 
  facet_wrap(~Env, scales = 'free_y')+
  ylab('Slope of trial regression') +
  theme_bw()+
  scale_x_discrete(labels = NULL) +
  theme(legend.position = 'bottom')
)


pdf(here(plot_dir, 'Optimum_posteriors.pdf'))
for(plot in optimum_posteriors) {
  print(plot)
}
print(slope_posterior_plot)
dev.off()

getTrialPosteriors = function(file) {
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  trial_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  trial_data = trial_data[order(trial_data$Env_trial),]
  ribbons = c()
  points = c()
  for(exp in trial_data$Experimento) {
    d_exp = subset(d,Experimento==exp)
    d_exp$y_resid = resid(lm(y~Tester,d_exp))
    X = seq(min(d_exp$Env_genotype_std),max(d_exp$Env_genotype_std),length=100)
    newdata=data.frame(
      Experimento=exp, ExpTester = NA,
      Env_genotype_std = X,
      Env_trial_std = d_exp$Env_trial_std[1]
    )
    # newdata$ExpTester[-1] = NA
    pre = posterior_epred(bm,newdata = newdata)
    pre_summary = apply(pre,2,function(y) {
      c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
    })
    pre_summary = data.frame(Env_genotype = X,t(pre_summary))
    rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
    pre_summary$Env_genotype = rescale(pre_summary$Env_genotype,bm_list)
    # pre_summary$low.2.5. = rescale(pre_summary$low.2.5.,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    ribbons = rbind(ribbons,data.frame(Experimento = exp,pre_summary))
    # points = rbind(points,data.frame(Experimento = exp,y = d_exp$y,Env_genotype = d_exp$Env_genotype))
    points = rbind(points,data.frame(Experimento = exp,y = d_exp$y_resid, Env_genotype = d_exp$Env_genotype))
  }
  ribbons$Experimento = factor(ribbons$Experimento,levels = trial_data$Experimento)
  ribbons_max = tapply(1:nrow(ribbons),ribbons$Experimento,function(x) ribbons$Env_genotype[x][which.max(ribbons$mean[x])])
  ribbons_max = data.frame(Experimento = factor(names(ribbons_max),levels = trial_data$Experimento),Env_genotype = ribbons_max)
  points$Experimento = factor(points$Experimento,levels = trial_data$Experimento)
  p=ggplot(ribbons) + 
    facet_wrap('Experimento',scales = 'free_y')+ 
    ggtitle(sprintf('%s, %s', trait, env)) + 
    xlab(env) + 
    ylab('scaled BLUP residual') +
    geom_ribbon(aes(x=Env_genotype, ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
    # geom_point(data = points,aes(x = Env_genotype,y=y)) + 
    geom_vline(data = ribbons_max,aes(xintercept = Env_genotype),color='red') + 
    geom_line(aes(x=Env_genotype,y=mean),linewidth=2,color='blue') +
    theme(text = element_text(size = 5)) +
    theme_bw()
  # print(p)
  p
}

experimento_posteriors = lapply(files, getTrialPosteriors)

pdf(here(plot_dir,'Experimento_posteriors.pdf'))
for(plot in experimento_posteriors) {
  print(plot)
}
dev.off()

getTrialPosteriorsCombined = function(file) {
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  trial_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  trial_data = trial_data[order(trial_data$Env_trial),]
  ribbons = c()
  points = c()
  for(exp in trial_data$Experimento) {
    d_exp = subset(d,Experimento==exp)
    d_exp$y_resid = resid(lm(y~Tester,d_exp))
    X = seq(min(d_exp$Env_genotype_std),max(d_exp$Env_genotype_std),length=100)
    newdata=data.frame(
      Experimento=exp, ExpTester = NA,
      Env_genotype_std = X,
      Env_trial_std = d_exp$Env_trial_std[1]
    )
    # newdata$ExpTester[-1] = NA
    pre = posterior_epred(bm,newdata = newdata)
    pre_summary = apply(pre,2,function(y) {
      c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
    })
    pre_summary = data.frame(Env_genotype = X,t(pre_summary))
    rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
    pre_summary$Env_genotype = rescale(pre_summary$Env_genotype,bm_list)
    # pre_summary$low.2.5. = rescale(pre_summary$low.2.5.,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    ribbons = rbind(ribbons,data.frame(Experimento = exp,pre_summary, Env_trial = d_exp$Env_trial[1]))
    # points = rbind(points,data.frame(Experimento = exp,y = d_exp$y,Env_genotype = d_exp$Env_genotype))
    points = rbind(points,data.frame(Experimento = exp,y = d_exp$y_resid, Env_genotype = d_exp$Env_genotype))
  }
  ribbons$Experimento = factor(ribbons$Experimento,levels = trial_data$Experimento)
  ribbons_max = tapply(1:nrow(ribbons),ribbons$Experimento,function(x) ribbons$Env_genotype[x][which.max(ribbons$mean[x])])
  ribbons_max = data.frame(Experimento = factor(names(ribbons_max),levels = trial_data$Experimento),Env_genotype = ribbons_max, Env_trial = trial_data$Env_trial)
  points$Experimento = factor(points$Experimento,levels = trial_data$Experimento)
  p=ggplot(ribbons) + 
    ggtitle(sprintf('%s, %s', trait, env)) + 
    xlab(sprintf('Accession %s-of-origin', env)) + 
    ylab('scaled BLUP residual') +
    # geom_ribbon(aes(x=Env_genotype,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
    # geom_point(data = points,aes(x = Env_genotype,y=y)) + 
    geom_vline(data = ribbons_max,aes(xintercept = Env_genotype, color = Env_trial)) + 
    geom_line(aes(x=Env_genotype,y=mean, group = Experimento,color = Env_trial),linewidth=1) +
    # scale_color_continuous(name = 'Trial value', low = 'blue', high = 'yellow') +
    scale_color_viridis(name = sprintf('trial %s', env)) +
    theme(text = element_text(size = 5)) +
    theme_bw()
  # print(p)
  p
}

experimento_posteriors_combined = lapply(files, getTrialPosteriorsCombined)
# experimento_posteriors_combined

pdf(here(plot_dir,'Experimento_posteriors_combined.pdf'))
for(plot in experimento_posteriors_combined) {
  print(plot)
}
dev.off()

alter_curve_plot = function(p) {
  env_var = p$labels$title
  units = ''
  if(grepl('elevation', env_var, fixed = T))
  {units = 'Trial elevation (m)'} else if(grepl('precipitation', env_var, fixed = T))
  {units = 'Trial precipitation (mm)'} else if(grepl('temperature', env_var, fixed = T))
  {units = 'Trial temperature (°C)'}
  p + 
    theme(legend.position = 'bottom',
          text = element_text(size = 6),
          legend.key.height = unit(.3, "lines"),
          legend.key.width = unit(1, "lines"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,0,0)) +
    scale_color_viridis(name = units) +
    ggtitle('')
}

alter_fit_plot = function(p) {
  p + 
    theme(legend.position = 'none',
          text = element_text(size = 6)) +
    ggtitle('')
}

(transfer_plots_brms = plot_grid(plot_grid(plotlist = lapply(experimento_posteriors_combined[c(2, 6, 10)], 
                                                             alter_curve_plot), 
                                                             ncol = 3, labels = c('A', 'B', 'C')),
                            plot_grid(plotlist = lapply(optimum_posteriors[c(2, 6, 10)], 
                                                        alter_fit_plot), 
                                      ncol = 3, labels = c('D', 'E', 'F')),
                            nrow = 2, rel_heights = c(1, 1)))

ggsave('transfer_plots_brms.pdf',
       plot = transfer_plots_brms,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 160,
       height = 80,
       units = 'mm'
)
       
ggsave('slope_posterior_plot.pdf',
       plot = slope_posterior_plot,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 3,
       height = 3,
       units = 'in'
)

(all_trial_posteriors = plot_grid(plotlist = lapply(experimento_posteriors_combined[c(1, 5, 9,
                                                                               2, 6, 10,
                                                                               3, 7, 11,
                                                                               4, 8, 12)],
                                                    function(p) p + theme(legend.position = 'bottom',
                                                                          text = element_text(size = 6),
                                                                          legend.key.height = unit(.3, "lines"),
                                                                          legend.key.width = unit(1, "lines"))), 
                                  ncol = 3, labels = 'auto'))

ggsave('all_trial_posteriors.pdf',
       plot = all_trial_posteriors,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 7,
       height = 7,
       units = 'in'
)

(all_optimum_posteriors = plot_grid(plotlist = lapply(optimum_posteriors[c(1, 5, 9,
                                                                            2, 6, 10,
                                                                            3, 7, 11,
                                                                            4, 8, 12)], 
                                                      function(p) p + theme(legend.position = 'none',
                                                                            text = element_text(size = 6),
                                                                            legend.key.height = unit(.3, "lines"),
                                                                            legend.key.width = unit(1, "lines"))), 
                                                      ncol = 3, labels = 'auto'))


ggsave('all_optimum_posteriors.pdf',
       plot = all_optimum_posteriors,
       device = 'pdf',
       here(plot_dir, 'Manuscript'),
       width = 7,
       height = 7,
       units = 'in'
)


## for talk
lapply(experimento_posteriors_combined[c(2, 6, 10)], 
       alter_curve_plot)
alter_curve_plot(experimento_posteriors_combined[[1]])

{
  file = files[2]
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  trial_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  trial_data = trial_data[order(trial_data$Env_trial),]
  ribbons = c()
  points = c()
  for(exp in trial_data$Experimento) {
    d_exp = subset(d,Experimento==exp)
    d_exp$y_resid = resid(lm(y~Tester,d_exp))
    X = seq(min(d_exp$Env_genotype_std),max(d_exp$Env_genotype_std),length=100)
    newdata=data.frame(
      Experimento=exp, ExpTester = NA,
      Env_genotype_std = X,
      Env_trial_std = d_exp$Env_trial_std[1]
    )
    # newdata$ExpTester[-1] = NA
    pre = posterior_epred(bm,newdata = newdata)
    pre_summary = apply(pre,2,function(y) {
      c(low = quantile(y,.025),high=quantile(y,.975),mean=mean(y))
    })
    pre_summary = data.frame(Env_genotype = X,t(pre_summary))
    rescale = function(x,bm_list) x*bm_list$std_x_sd + bm_list$std_x_mean
    pre_summary$Env_genotype = rescale(pre_summary$Env_genotype,bm_list)
    # pre_summary$low.2.5. = rescale(pre_summary$low.2.5.,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    # pre_summary$geno_x = rescale(pre_summary$geno_x,bm_list)
    ribbons = rbind(ribbons,data.frame(Experimento = exp,pre_summary, Env_trial = d_exp$Env_trial[1]))
    # points = rbind(points,data.frame(Experimento = exp,y = d_exp$y,Env_genotype = d_exp$Env_genotype))
    points = rbind(points,data.frame(Experimento = exp,y = d_exp$y_resid, Env_genotype = d_exp$Env_genotype))
  }
  ribbons$Experimento = factor(ribbons$Experimento,levels = trial_data$Experimento)
  ribbons_max = tapply(1:nrow(ribbons),ribbons$Experimento,function(x) ribbons$Env_genotype[x][which.max(ribbons$mean[x])])
  ribbons_max = data.frame(Experimento = factor(names(ribbons_max),levels = trial_data$Experimento),Env_genotype = ribbons_max, Env_trial = trial_data$Env_trial)
  points$Experimento = factor(points$Experimento,levels = trial_data$Experimento)
  p=ggplot(ribbons %>% filter(Experimento == 'mBCL2011')) + 
    # ggtitle(sprintf('%s, %s', trait, env)) + 
    xlab(sprintf('Accession %s-of-origin', env)) + 
    ylab('Accession field weight (scaled BLUP)') +
    # geom_ribbon(aes(x=Env_genotype,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
    geom_point(data = points %>% filter(Experimento == 'mBCL2011'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max ,aes(xintercept = Env_genotype, color = Env_trial)) + 
    geom_line(aes(x=Env_genotype,y=mean, group = Experimento,color = Env_trial),linewidth=1) +
    # scale_color_continuous(name = 'Trial value', low = 'blue', high = 'yellow') +
    scale_color_viridis(name = 'Trial elevation (m)') +
    theme(text = element_text(size = 5)) +
    theme_bw()
  # print(p)
  p
}

## for poster
# (all_trial_posteriors = plot_grid(plotlist = lapply(experimento_posteriors_combined[c(2, 6, 10,
#                                                                                       3, 7, 11)],
#                                                     function(p) p + theme(legend.position = 'bottom') + ggtitle('')), 
#                                   ncol = 3, labels = 'auto'))


####################################################################################
## manuscript plotting


# (environmental_map_transfer_plots = plot_grid(
#   plot_grid(elevation_map, transfer_plots, labels = 'AUTO', nrow = 2, rel_heights = c(2,1))))
# png(here(plot_dir, 'Manuscript', 'environmental_map_transfer_plots.png'), width = 600, height = 600, dpi = 300)
# print(environmental_map_transfer_plots)
# dev.off()

## old manual plotting
# for(trait in traitNames) {
#   # print(data %>% filter(SampleID %in% cimmyt_grow$Sample.ID) %>% pull(SampleID) %>% unique() %>% length())
#   ## plot transfer distance plots
#   data_trait = subset(data %>% filter(SampleID %in% cimmyt_grow$Sample.ID),Trait==trait & !is.na(Value) & !is.na(Elevation_distance))
#   data_trait$y = resid(lm(Value~Experimento:Tester,data_trait)) * 100
#   elevation_plots[[trait]] = ggplot(data_trait,aes(x=x_axis,y = y, color = Elevation_trial)) +
#     geom_smooth(aes(group = Experimento,color = Elevation_trial),se=F,method = 'lm',formula = plot_formula) +
#     scale_color_continuous(name = 'Trial elevation', low = 'blue', high = 'red')+
#     ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
#     theme_bw() +
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 10)) +
#     xlab('Transfer distance') +
#     ylab('Trait residual')
#   temp_plots[[trait]] = ggplot(data_trait,aes(x=x_axis,y = y, color = meanTemp_trial)) +
#     geom_smooth(aes(group = Experimento,color = meanTemp_trial),se=F,method = 'lm',formula = plot_formula) +
#     scale_color_continuous(name = 'Trial temperature', low = 'blue', high = 'red')+
#     ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
#     theme_bw() +
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 10)) +
#     xlab('Transfer distance') +
#     ylab('Trait residual')
#   precipitation_plots[[trait]] = ggplot(data_trait,aes(x=abs(annualPrecipitation_distance),y = y, color = annualPrecipitation_trial)) +
#     geom_smooth(aes(group = Experimento,color = annualPrecipitation_trial),se=F,method = 'lm',formula = plot_formula) +
#     scale_color_continuous(name = 'Trial precipitation', low = 'blue', high = 'red')+
#     ggtitle(gsub("([[:upper:]])", " \\1", trait), ) +
#     theme_bw() +
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 10)) +
#     xlab('Transfer distance') +
#     ylab('Trait residual')
#   
#   ## modeling with average per trial, or using genotype as random effect - did not see sig difference in optimal
#   # m = lmer(y~poly(Env_genotype,2):Experimento + Experimento:Tester + (1|SampleID),data_trait)
#   # m = lm(y~poly(Env_genotype,2):Experimento + Experimento:Tester,data_trait)
#   # coefs_m = fixef(m)
#   # anova_m = summary(m)
#   # 
#   # coefs = do.call(rbind,lapply(unique(data_trait$Experimento),function(exp) {
#   #   i = data_trait$Experimento == exp
#   #   t1 = data_trait$Tester[i][1]
#   #   ID1 = data_trait$SampleID[i][1]
#   #   x = do.call(seq,as.list(c(range(data_trait[i,]$Env_genotype),length=1000)))
#   #   predictions_trial = predict(m,newdata = data.frame(Env_genotype=x,Tester=t1,Experimento = exp,SampleID=ID1))
#   #   data.frame(Experimento = exp, 
#   #              X2 = coefs_m[sprintf('poly(Env_genotype, 2)1:Experimento%s',exp)],
#   #              p = anova_m$coefficients[sprintf('poly(Env_genotype, 2)1:Experimento%s',exp),5],
#   #              Max_X = x[which.max(predictions_trial)])
#   # }))
#   
#   ## calculate coefficients of parabolas and create models
#   elevation_coefs = as.data.frame(do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) {
#     t1 = data_trait$Tester[i][1]
#     m = lm(y~poly(Elevation_genotype,2)+Tester,data_trait[i,])
#     # m2 = lm(y~poly(Env_genotype,2):Tester+Latitude_genotype+Tester,data_trait[i,])#
#     x = do.call(seq,as.list(c(range(data_trait[i,]$Elevation_genotype),length=1000)))
#     if(coef(m)[3] > 0) {
#       c(coef(m)[1:3],x[which.min(predict(m,newdata = data.frame(Elevation_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     } else {
#       c(coef(m)[1:3],x[which.max(predict(m,newdata = data.frame(Elevation_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     }
#   })))
#   colnames(elevation_coefs)[1:4] = c('Intercept','X','X2','vertex')
#   elevation_models[[trait]] = data.frame(
#     Trait = trait,
#     Env = tapply(data_trait$Elevation_trial,data_trait$Experimento,mean),
#     elevation_coefs
#     # do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) coef(lm(y~poly(Elevation_distance,2),data_trait[i,]))))
#   )
#   
#   temp_coefs = as.data.frame(do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) {
#     t1 = data_trait$Tester[i][1]
#     m = lm(y~poly(meanTemp_genotype,2)+Tester,data_trait[i,])
#     # m2 = lm(y~poly(Env_genotype,2):Tester+Latitude_genotype+Tester,data_trait[i,])#
#     x = do.call(seq,as.list(c(range(data_trait[i,]$meanTemp_genotype),length=1000)))
#     if(coef(m)[3] > 0) {
#       c(coef(m)[1:3],x[which.min(predict(m,newdata = data.frame(meanTemp_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     } else {
#       c(coef(m)[1:3],x[which.max(predict(m,newdata = data.frame(meanTemp_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     }
#   })))
#   colnames(temp_coefs)[1:4] = c('Intercept','X','X2','vertex')
#   temp_models[[trait]] = data.frame(
#     Trait = trait,
#     Env = tapply(data_trait$meanTemp_trial,data_trait$Experimento,mean),
#     temp_coefs
#     # do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) coef(lm(y~poly(meanTemp_distance,2),data_trait[i,]))))
#   )
#   
#   prec_coefs = as.data.frame(do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) {
#     t1 = data_trait$Tester[i][1]
#     m = lm(y~poly(annualPrecipitation_genotype,2)+Tester,data_trait[i,])
#     # m2 = lm(y~poly(Env_genotype,2):Tester+Latitude_genotype+Tester,data_trait[i,])#
#     x = do.call(seq,as.list(c(range(data_trait[i,]$annualPrecipitation_genotype),length=1000)))
#     if(coef(m)[3] > 0) {
#       c(coef(m)[1:3],x[which.min(predict(m,newdata = data.frame(annualPrecipitation_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     } else {
#       c(coef(m)[1:3],x[which.max(predict(m,newdata = data.frame(annualPrecipitation_genotype=x,Tester=t1)))],unlist(anova(m)[1,4:5]))
#     }
#   })))
#   colnames(prec_coefs)[1:4] = c('Intercept','X','X2','vertex')
#   precipitation_models[[trait]] = data.frame(
#     Trait = trait,
#     Env = tapply(data_trait$annualPrecipitation_trial,data_trait$Experimento,mean),
#     prec_coefs
#     # do.call(rbind,tapply(1:nrow(data_trait),data_trait$Experimento,function(i) coef(lm(y~poly(annualPrecipitation_distance,2),data_trait[i,]))))
#   )
# }

# ggplot(models,aes(x=Env,y=X)) + geom_point() + facet_wrap(~Trait,scales = 'free')
# ggplot(models,aes(x=Env/1000,y=X2)) + geom_point() + facet_wrap(~Trait,scales = 'free') + xlab('Elevation (100m) of trial') + ylab('Coefficient of (sample-trial elevation)^2')
# 
# 
# plots2 = list()
# models2 = list()
# for(trait in traitNames) {
#   data_trait = subset(data,Trait==trait & !is.na(Value) & !is.na(Env_distance))
#   # data_trait$y = resid(lm(Value~Experimento:Tester,data_trait))
#   data_trait$y = data_trait$Value
#   data_trait$y_resid = resid(lm(Value~Experimento:Tester,data_trait))
# 
#   plots2[[trait]] =
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
# (p = cowplot::plot_grid(plotlist = plots2,ncol=2))
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
