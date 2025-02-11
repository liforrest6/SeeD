# script to generate brms plots for transfer distance

library(cowplot)
library(data.table)
library(MASS)
library(brms)
library(viridis)
library(scales)

source(here::here('config.R'))

##### generate brms models for transfer plots to model ####################################################################################
## do not use if generating models
## from Optimal_model_environment_unconstrained.R on farm


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

## read brms model data here
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

#### plotting models#######################################################################

## read constrained files
files = list.files(here(analyses_dir, 'TransferPlots'), pattern='.rds',full.names = T)
## read unconstrained files
files = list.files(here(analyses_dir, 'LocalAdaptation/Results_unconstrained'), pattern='.rds',full.names = T)

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

## barplot of b_h for each trial for three environmental variables, and four yield traits
(slope_posterior_plot = ggplot(slope_posteriors,aes(x=Trait,y=b_h_x_trial, color = Trait)) + 
    geom_boxplot(aes(group = Trait)) + 
    geom_hline(yintercept=0) + 
    facet_wrap(~Env, scales = 'free_y')+
    ylab('Slope of trial regression') +
    theme_bw()+
    scale_x_discrete(labels = NULL) +
    theme(legend.position = 'bottom')
)

## see local adaptation curve for each trial
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
    geom_point(data = points, aes(x = Env_genotype, y = y), alpha = 0.2, size = 0.25) +
    geom_vline(data = ribbons_max,aes(xintercept = Env_genotype),color='red') + 
    geom_line(aes(x=Env_genotype,y=mean),linewidth=2,color='blue') +
    theme(text = element_text(size = 5)) +
    theme_bw()
  # print(p)
  p
}

## generate transfer plot curves for each trial so we can put them in large pdf
experimento_posteriors = lapply(files, getTrialPosteriors)

## get transfer plot curves for each trial, but combine them for env:trait for combined viewing
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

## formatting function
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

## formatting function
alter_fit_plot = function(p) {
  p + 
    theme(legend.position = 'none',
          text = element_text(size = 6)) +
    ggtitle('')
}

optimum_posteriors$`/Users/liforrest/Documents/Projects/SeeD/Analyses//TransferPlots/bm_precipitation_FieldWeight.rds` + ylim(-100, NA)

## generates manuscript figure with both combined posteriors and regression
(transfer_plots_brms = plot_grid(plot_grid(plotlist = lapply(experimento_posteriors_combined[c(2, 6, 10)], 
                                                             alter_curve_plot), 
                                           ncol = 3, labels = c('A', 'B', 'C')),
                                 plot_grid(plotlist = lapply(optimum_posteriors[c(2, 6, 10)], 
                                                             alter_fit_plot), 
                                           ncol = 3, labels = c('D', 'E', 'F')),
                                 nrow = 2, rel_heights = c(1, 1)))


## generates figure looking at all env:trait combinations with all trials in each panel
(all_trial_posteriors = plot_grid(plotlist = lapply(experimento_posteriors_combined[c(1, 5, 9,
                                                                                      2, 6, 10,
                                                                                      3, 7, 11,
                                                                                      4, 8, 12)],
                                                    function(p) p + theme(legend.position = 'bottom',
                                                                          text = element_text(size = 6),
                                                                          legend.key.height = unit(.3, "lines"),
                                                                          legend.key.width = unit(1, "lines"))), 
                                  ncol = 3, labels = 'auto'))


## generates figure looking at all env:trait combinations with regression of trials in each panel
(all_optimum_posteriors = plot_grid(plotlist = lapply(optimum_posteriors[c(1, 5, 9,
                                                                           2, 6, 10,
                                                                           3, 7, 11,
                                                                           4, 8, 12)], 
                                                      function(p) p + theme(legend.position = 'none',
                                                                            text = element_text(size = 6),
                                                                            legend.key.height = unit(.3, "lines"),
                                                                            legend.key.width = unit(1, "lines"))), 
                                    ncol = 3, labels = 'auto'))


## generates figure with regression of optimum posteriors for all trials, separated by env:trait
pdf(here(plot_dir, 'Optimum_posteriors.pdf'))
for(plot in optimum_posteriors) {
  print(plot)
}
print(slope_posterior_plot)
dev.off()

## generates figure for posteriors for every trial, separately (the big one)
pdf(here(plot_dir, 'Experimento_posteriors.pdf'))
for(plot in experimento_posteriors) {
  print(plot)
}
dev.off()

## generates figure by putting trials together for env:trait combination
pdf(here(plot_dir, 'Experimento_posteriors_combined.pdf'))
for(plot in experimento_posteriors_combined) {
  print(plot)
}
dev.off()

## generates manuscript figure
ggsave('transfer_plots_brms.tif',
       plot = transfer_plots_brms,
       device = 'tiff',
       here(plot_dir, 'Manuscript'), ## change for unconstrained
       width = 160,
       height = 80,
       units = 'mm'
)

## generates regression slope
ggsave('slope_posterior_plot.pdf',
       plot = slope_posterior_plot,
       device = 'pdf',
       here(plot_dir, 'Manuscript'), ## change pdf location
       width = 3,
       height = 3,
       units = 'in'
)

## manuscript transfer distance plots for all trials across all traits
ggsave('all_trial_posteriors.pdf',
       plot = all_trial_posteriors,
       device = 'pdf',
       here(plot_dir, 'Manuscript'), ## change location of pdf
       width = 7,
       height = 7,
       units = 'in'
)

## all optimums across all traits
ggsave('all_optimum_posteriors.pdf',
       plot = all_optimum_posteriors,
       device = 'pdf',
       here(plot_dir, 'unconstrained'), ## change location of pdf
       width = 7,
       height = 7,
       units = 'in'
)

## read independent files
files = list.files(here(analyses_dir, 'LocalAdaptation/Results_independent'), pattern='.rds',full.names = T)
independent_optimum_posteriors = list()
independent_slope_posteriors = c()
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
  h_posteriors_tall$geno_x = bm_list$data_trait$geno_x_std[match(h_posteriors_tall$name,bm_list$data_trait$SampleID)]
  
  lines = as_draws_df(bm,c('b_h_Intercept'))
  # slope_posteriors = rbind(slope_posteriors,data.frame(Trait = trait,b_h_x_trial = lines$b_h_Env_trial_std))
  independent_slope_posteriors = rbind(independent_slope_posteriors,data.frame(Trait = trait,Env=env))
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
  independent_optimum_posteriors[[file]] = p
}

getTrialPosteriors_independent = function(file) {
  trait = strsplit(file,'_')[[1]][[4]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[3]]
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
    ## generate points for outcome variability
    geom_point(data = points, aes(x = Env_genotype, y = y), alpha = 0.2, size = 0.25) +
    # geom_point(data = points,aes(x = Env_genotype,y=y)) + 
    geom_vline(data = ribbons_max,aes(xintercept = Env_genotype),color='red') + 
    geom_line(aes(x=Env_genotype,y=mean),linewidth=2,color='blue') +
    theme(text = element_text(size = 5)) +
    theme_bw()
  # print(p)
  p
}

independent_experimento_posteriors = lapply(files, getTrialPosteriors_independent)

## generates figure for posteriors for every trial, separately (the big one) for unconstrained 
pdf(here(plot_dir, 'unconstrained', 'Experimento_posteriors.pdf'))
for(plot in independent_experimento_posteriors) {
  print(plot)
}
dev.off()




#### generates single trial plots for manuscript walk-through #### 
# lapply(experimento_posteriors_combined[c(2, 6, 10)],
#        alter_curve_plot)
# alter_curve_plot(experimento_posteriors_combined[[1]])

{
  file = files[2]
  trait = strsplit(file,'_')[[1]];trait = trait[length(trait)];trait = sub('.rds','',trait,fixed=T)
  env = strsplit(file,'_')[[1]][[2]]
  bm_list = readRDS(file)
  bm = bm_list$bm
  d = bm_list$data_trait
  trial_data = droplevels(d[sapply(unique(d$Experimento),function(x) which(d$Experimento == x)[1]),])
  trial_data = trial_data[order(trial_data$Env_trial),]
  
  ## for generating slope posterior plot
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
  # h_posteriors_summary$optimum_genotype = merge(h_posteriors_summary[c('name')], ribbons_max[c('Experimento', 'Env_genotype')], by.x = 'name', by.y = 'Experimento')$Env_genotype
  
  
  ## for generating trials
  ribbons = c()
  points = c()
  posterior_draws = c()
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
    
    ### to draw posterior draws
    pre_draws = t(pre[1:1000,])
    pre_draws_df = data.frame(Env_genotype = X, pre_draws)
    pre_draws_df$Env_genotype = rescale(pre_draws_df$Env_genotype, bm_list)
    posterior_draws = rbind(posterior_draws, data.frame(Experimento = exp, Env_trial = d_exp$Env_trial[1], pre_draws_df))
    
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
  # h_posteriors_summary$optimum_genotype = merge(h_posteriors_summary[c('name')], ribbons_max[c('Experimento', 'Env_genotype')], by.x = 'name', by.y = 'Experimento')$Env_genotype
  ribbons_max$apex_posterior_mean = h_posteriors_summary[match(ribbons_max$Experimento, h_posteriors_summary$name), 'mean']
  
  points$Experimento = factor(points$Experimento,levels = trial_data$Experimento)
  
  posterior_draws$Experimento = factor(posterior_draws$Experimento, levels = trial_data$Experimento)
  posterior_draws_tall = pivot_longer(posterior_draws, cols = starts_with('X'), names_to = 'draw', values_to = 'value')
  

}


(
  m2012BCL_transfer=ggplot(ribbons %>% filter(Experimento == 'm2012BCL')) +
  # ggtitle(sprintf('%s, %s', trait, env)) +
  xlab(sprintf('Accession %s-of-origin', env)) +
  ylab('Accession field weight (scaled BLUP)') +
  geom_ribbon(aes(x=Env_genotype,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
  geom_point(data = points %>% filter(Experimento == 'm2012BCL'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012BCL'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
    
    # geom_line(aes(x=Env_genotype,y=mean, group = Experimento,color = Env_trial),linewidth=1) +
    scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#619CFF')) +
    theme_bw()+ 
    theme(legend.position = 'bottom') +
    ggtitle('m2012BCL')
)

(
  m2012AAFTU_transfer=ggplot(ribbons %>% filter(Experimento == 'm2012AAFTU')) +
  # ggtitle(sprintf('%s, %s', trait, env)) +
  xlab(sprintf('Accession %s-of-origin', env)) +
  ylab('Accession field weight (scaled BLUP)') +
  geom_ribbon(aes(x=Env_genotype,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
  geom_point(data = points %>% filter(Experimento == 'm2012AAFTU'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012AAFTU'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
    
    # geom_line(aes(x=Env_genotype,y=mean, group = Experimento,color = Env_trial),linewidth=1) +
  scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#F564E3')) +
    theme_bw()+ 
    theme(legend.position = 'bottom') +
    ggtitle('m2012AAFTU')
)

(
  m2012BNY_transfer=ggplot(ribbons %>% filter(Experimento == 'm2012BNY')) +
    # ggtitle(sprintf('%s, %s', trait, env)) +
    xlab(sprintf('Accession %s-of-origin', env)) +
    ylab('Accession field weight (scaled BLUP)') +
    geom_ribbon(aes(x=Env_genotype,ymin = low.2.5.,ymax = high.97.5.),alpha = 0.5) +
    geom_point(data = points %>% filter(Experimento == 'm2012BNY'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012BNY'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
    
    # geom_line(aes(x=Env_genotype,y=mean, group = Experimento,color = Env_trial),linewidth=1) +
    scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#00BA38')) +
    theme_bw()+ 
    theme(legend.position = 'none') +
    ggtitle('m2012BNY')
)


(
  m2012BCL_posterior=ggplot(posterior_draws_tall %>% filter(Experimento == 'm2012BCL')) +
    # ggtitle(sprintf('%s, %s', trait, env)) +
    xlab(sprintf('Accession %s-of-origin', env)) +
    ylab('Accession field weight (scaled BLUP)') +

    geom_line(aes(x=Env_genotype,y=value, group = draw),linewidth=1, alpha = 0.05) +

    geom_point(data = points %>% filter(Experimento == 'm2012BCL'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012BCL'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
    scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#619CFF')) +
    theme_bw()+ 
    theme(legend.position = 'none') +
    ggtitle('m2012BCL')
)

(
  m2012BNY_posterior=ggplot(posterior_draws_tall %>% filter(Experimento == 'm2012BNY')) +
    # ggtitle(sprintf('%s, %s', trait, env)) +
    xlab(sprintf('Accession %s-of-origin', env)) +
    ylab('Accession field weight (scaled BLUP)') +

    geom_line(aes(x=Env_genotype,y=value, group = draw),linewidth=1, alpha = 0.05) +

    geom_point(data = points %>% filter(Experimento == 'm2012BNY'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012BNY'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
    scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#F564E3')) +
    theme_bw()+ 
    theme(legend.position = 'none') +
    ggtitle('m2012BNY')
)


(
  m2012AAFTU_posterior=ggplot(posterior_draws_tall %>% filter(Experimento == 'm2012AAFTU')) +
    # ggtitle(sprintf('%s, %s', trait, env)) +
    xlab(sprintf('Accession %s-of-origin', env)) +
    ylab('Accession field weight (scaled BLUP)') +

    geom_line(aes(x=Env_genotype,y=value, group = draw),linewidth=1, alpha = 0.05) +
    

    geom_point(data = points %>% filter(Experimento == 'm2012AAFTU'),aes(x = Env_genotype,y=y)) +
    geom_vline(data = ribbons_max %>% filter(Experimento == 'm2012AAFTU'),
               aes(xintercept = apex_posterior_mean, color = factor(apex_posterior_mean))) +
        scale_color_manual(name = 'Optimum collection elevation (masl)', values = c('#00BA38')) +
    theme_bw()+ 
    theme(legend.position = 'none') +
    ggtitle('m2012AAFTU')
)

(
annotated_optimum = ggplot() + ggtitle(sprintf('%s, %s', trait, env)) + 
                       xlab(sprintf('Trial %s', env)) +
                       ylab(sprintf('Optimum %s', env)) +
                       geom_ribbon(data = fit_df,aes(x=x,ymin = low,ymax = high),alpha = 0.3)+
                       geom_line(data = fit_df,aes(x=x,y=mean))+
                       geom_abline(slope=1,intercept=0,color='red')+
                       geom_segment(data=h_posteriors_summary,aes(x=Env_trial,y=`low.2.5.`,xend=Env_trial,yend=`high.97.5.`,group=name,color=name))+
                       geom_point(data=h_posteriors_summary,aes(x=Env_trial,y=mean,color=name),size=2) +
                       theme_bw() +
    geom_hline(data = h_posteriors_summary %>% 
                   filter(name %in% (c('m2012AAFTU',
                                        'm2012BNY', 
                                        'm2012BCL'))),
                 aes(yintercept = mean, color = name)) +
    geom_text(aes(0,filter(h_posteriors_summary, name == 'm2012AAFTU')$mean,label = 'A', vjust = -1)) +
    geom_text(aes(0,filter(h_posteriors_summary, name == 'm2012BNY')$mean,label = 'B', vjust = -1)) +
    geom_text(aes(0,filter(h_posteriors_summary, name == 'm2012BCL')$mean,label = 'C', vjust = -1)) 
    
                     # geom_jitter(alpha = 0.01,aes(y=value,color=name)) +
                     # geom_violin(aes(y=value,group = name,color=name),position = position_dodge())
)


(
  posterior_draw_combined_plot = plot_grid(plot_grid(m2012AAFTU_posterior, m2012BNY_posterior, m2012BCL_posterior, nrow = 1, labels = c('A', 'B', 'C')),
                                  annotated_optimum, labels = c('', 'D'), nrow = 2)
)

(
  transfer_optimum_combined_plot = plot_grid(plot_grid(m2012AAFTU_transfer, m2012BNY_transfer, m2012BCL_transfer, nrow = 1, labels = c('A', 'B', 'C')),
                                             annotated_optimum, labels = c('', 'D'), nrow = 2)
)

ggsave('posterior_draw_combined_plot.tif',
       plot = posterior_draw_combined_plot,
       device = 'tiff',
       here(plot_dir, 'Manuscript'), ## change pdf location
       width = 8,
       height = 8,
       units = 'in'
)
