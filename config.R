####################################################################################
# Configuration file for SeeDs GEA analysis
#
# Author: Forrest Li
# Load all libraries, source scripts
####################################################################################

library(readxl)
library(dplyr)
library(stringr)
library(data.table)
library(here)
library(tidyr)
library(purrr)



here::here()
genetic_data_dir = paste0(here::here(), '/Genetic_data/')
phenotype_data_dir = paste0(here::here(), '/Phenotype_data/')
env_data_dir = paste0(here::here(), '/Env_data/')

plot_dir = paste0(here::here(), '/Plots/')
scripts_dir = paste0(here::here(), '/Scripts/')
analyses_dir = paste0(here::here(), '/Analyses/')

source(here::here("Scripts/analysisFunctions.R"))
source(here::here("Scripts/plotFunctions.R"))

