library(brms)
library(data.table)
setwd('/group/runciegrp2/Projects/SeeD/Analyses/LocalAdaptation')
outdir = 'Results_unconstrained'
try(dir.create(outdir,showWarnings = FALSE))

data = fread('../../Phenotype_data/blups_deregressed_CV.csv',data.table = F)
trial_info = read.csv('../../Phenotype_data/Trial_info.csv')
trial_worldclim = read.csv('../../Env_data/Trial_worldclim.csv')
trial_info = merge(trial_info, trial_worldclim, by = 'Experimento')

dir = '/group/runciegrp2/Projects/SeeD/'
elevation_df = fread(paste0(dir, 'Env_data/GEA-climate-nontransformed.csv'), data.table = F)
env_data = fread(paste0(dir, 'Env_data/SEEDGWAS_worldclim.csv'), data.table = F)
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

env_data = merge(env_data, elevation_df[, c('Unique.ID', 'elevation')], by = 'Unique.ID')

geno_info = merge(geno_info,env_data,by.x = 'V1',by.y = 'Unique.ID')


# geno_env_variable = 'elevation'
# trial_env_variable = 'Trial_elevation'
# trait = 'FieldWeight'

geno_env_variable = c('elevation','annualPrecipitation', 'meanTemp')[as.numeric(commandArgs(t=T)[1])]
trial_env_variable = c('Trial_elevation','annualPrecipitation', 'meanTemp')[as.numeric(commandArgs(t=T)[1])]
traits = c('FieldWeight', 'GrainWeightPerHectareCorrected', 'BareCobWeight', 'PlantHeight')
trait = traits[as.numeric(commandArgs(t=T)[2])]


data_trait = subset(data,Trait==trait)
data_trait$Env_genotype = geno_info[[geno_env_variable]][match(data_trait$SampleID,geno_info$Sample)]
data_trait$Env_trial = trial_info[[trial_env_variable]][match(data_trait$Experimento,trial_info$Experimento)]
data_trait = subset(data_trait,!is.na(Env_genotype))

std_x_mean = mean(data_trait$Env_genotype)
std_x_sd = sd(data_trait$Env_genotype)

data_trait$Env_genotype_std = (data_trait$Env_genotype - std_x_mean)/std_x_sd
data_trait$Env_trial_std = (data_trait$Env_trial - std_x_mean)/std_x_sd
data_trait$y = scale(data_trait$Value)

data_trait$ExpTester = interaction(data_trait$Experimento,data_trait$Tester,drop=TRUE)

## set priors, take out ub = 0 for parameter a to unconstrain from negative curve
prior1 = prior(normal(0,1),nlpar="c") + prior(normal(0,1),nlpar="a") + prior(normal(0,3),nlpar="h")
bm = brm(bf(y ~ c+a*(Env_genotype_std-h)^2, a~ 0+Experimento,h~1+Env_trial_std+(1|Experimento),sigma~0+Experimento,c~1+(1|ExpTester), nl = TRUE),
         data = data_trait, prior = prior1,
         control = list(adapt_delta = 0.9),cores=4)

saveRDS(list(bm=bm,data_trait=data_trait,std_x_mean=std_x_mean,std_x_sd=std_x_sd),file = file.path(outdir,sprintf('bm_%s_%s.rds',geno_env_variable,trait)))


