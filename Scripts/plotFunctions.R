####################################################################################
# Plotting functions for SeeDs GEA analysis
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################

library(ggplot2)
library(qqman)
library(ggmap)
library(ggpubr)
library(png)
library(grid)
library(gridExtra)

tmin_results = read_GEMMA_results(here(analyses_dir, 'GEA_output', 'GEMMA_univariate'), 'tmin', 'p_wald')
inv4mCoords = c(1.71e+8, 1.86e8)
hsftf9Coords = c(1.480e8 , 1.490e8)
inv4mSNPs = tmin_results %>% filter(BP > inv4mCoords[1] & BP < inv4mCoords[2] & Chr == 4) %>% pull(X)
hsftf9SNPs = tmin_results %>% filter(BP > hsftf9Coords[1] & BP < hsftf9Coords[2] & Chr == 9) %>% pull(X)