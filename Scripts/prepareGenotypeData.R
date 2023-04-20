##########################################
# Prepare genotype data
#
# Author: Forrest Li
# Functions for running analyses
##########################################

source(here::here('config.R'))

## 4692 SEEDGWAS with no repeat sequencing
non_dup_acc = read.csv(here(genetic_data_dir, 'selected_genotypeIDs.csv'))
## CIMMYT climate data
cimmyt = read_excel(here(env_data_dir, 'cimmyt_data.xlsx'), 
                    sheet = 'climate growing and fl season',
                    guess_max = 12000)
## Phenotyped genotypes mapping sheet
blups = read_excel(here(phenotype_data_dir, 'Blups_01.xlsx'), sheet = 'DNA to TC GID')

## 4020 unique Accession GID
blups$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)` %>% unique() %>% length()
## 4704 unique testcrosses
blups$`TC GID` %>% unique() %>% length()
## 4710 unique SEEDGWAS IDs
blups$`Sample ID`%>% unique() %>% length()

## 4692 SEEDGWAS genotypes in blups file
(non_dup_acc_GWAS_ID = blups[which(blups$`Sample ID`%in% non_dup_acc$Sample),])

## confirmed 18 with repeat sequencing
# blups[which(!blups$`Sample ID` %in% non_dup_acc$Sample),] 

## but only see one line for Sample ID, Accession ID, or TC GID
# blups[which(blups$`Sample ID` == 'SEEDGWAS3234'),] 
# blups[which(blups$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)`==24841),]
# blups[which(blups$`TC GID` == '238738'),]

## 2895 accessions in CIMMYT climate with GID matching the Accession GID from non-duplicate genotypes
finalAccessions = cimmyt[which(cimmyt$GID %in% non_dup_acc_GWAS_ID$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)`),]

## 3520 genotypes with climate data from CIMMYT
finalGenotypes = blups[which(blups$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)` %in% finalAccessions$GID),]
finalGenotypes

finalMat = merge(finalGenotypes[, 1:5], 
                 finalAccessions, 
                 by.x = 'AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)',
                 by.y = 'GID')


if(!file.exists(here(env_data_dir, 'GEA_finalGenotypeList.csv'))) {
  write.csv(finalMat, here(env_data_dir, 'GEA_finalGenotypeList.csv'))
}


