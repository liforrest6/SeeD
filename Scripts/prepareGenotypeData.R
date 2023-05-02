####################################################################################
# Process genotype set for GEA analyses
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################

source(here::here('config.R'))

## 4692 SEEDGWAS with no repeat sequencing
non_dup_acc = read.csv(here(genetic_data_dir, 'selected_genotypeIDs.csv'))
colnames(non_dup_acc) = c('Unique.ID', 'Sample')
## CIMMYT climate data
cimmyt = read_excel(here(env_data_dir, 'cimmyt_data.xlsx'), 
                    sheet = 'climate growing and fl season',
                    guess_max = 12000)
## Phenotyped genotypes mapping sheet
dna_to_tc_gid = read_excel(here(phenotype_data_dir, 'Blups_01.xlsx'), sheet = 'DNA to TC GID')
# blups = as.data.frame(readxl::read_xlsx(here(phenotype_data_dir, 'Blups_01.xlsx'),sheet = 'Blups01',col_types = c(rep('text',11),'numeric')))

## 4020 unique Accession GID
dna_to_tc_gid$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)` %>% unique() %>% length()
## 4704 unique testcrosses
dna_to_tc_gid$`TC GID` %>% unique() %>% length()
## 4710 unique SEEDGWAS IDs
dna_to_tc_gid$`Sample ID`%>% unique() %>% length()

## 4692 SEEDGWAS genotypes in blups file
(non_dup_acc_GWAS_ID = dna_to_tc_gid[which(dna_to_tc_gid$`Sample ID`%in% non_dup_acc$Sample),])
## there are five genotypes with no testcrosses whatsoever, keep for GEA anyways
bad_testcrosses = non_dup_acc_GWAS_ID[which(!non_dup_acc_GWAS_ID$`Sample ID` %in% (subset(non_dup_acc_GWAS_ID, `TC GID` != 286910) %>% pull(`Sample ID`))),]
## filter out TC GID 286910 for duplicate testcross, 4690
non_dup_acc_GWAS_ID = subset(non_dup_acc_GWAS_ID, `TC GID` != 286910 | is.na(`TC GID`))

## confirmed 18 with repeat sequencing
# dna_to_tc_gid[which(!dna_to_tc_gid$`Sample ID` %in% non_dup_acc$Sample),] 

## but only see one line for Sample ID, Accession ID, or TC GID
# dna_to_tc_gid[which(dna_to_tc_gid$`Sample ID` == 'SEEDGWAS3234'),] 
# dna_to_tc_gid[which(dna_to_tc_gid$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)`==24841),]
# dna_to_tc_gid[which(dna_to_tc_gid$`TC GID` == '238738'),]

## 2895 accessions in CIMMYT climate with GID matching the Accession GID from non-duplicate genotypes
accessions = cimmyt[which(cimmyt$GID %in% non_dup_acc_GWAS_ID$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)`),]

## 3516 genotypes with climate data from CIMMYT
genotypes = non_dup_acc_GWAS_ID[which(non_dup_acc_GWAS_ID$`AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)` %in% accessions$GID),]

genotypeList = merge(genotypes[, 1:5], 
                 accessions[, 1:6], 
                 by.x = 'AccessionGID (Germplasm ID at CIMMYT; General_identifier in Germinate)',
                 by.y = 'GID')

genotypeList = merge(genotypeList,
                     non_dup_acc,
                     by.x = 'Sample ID',
                     by.y = 'Sample')


if(file.exists(here(env_data_dir, 'GenotypesWithIDs.csv'))) {
  write.csv(genotypeList, here(env_data_dir, 'GenotypesWithIDs.csv'), row.names = F)
}


