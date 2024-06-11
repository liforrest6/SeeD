import numpy as np
import pandas as pd
import glob
import os

print('Start clumping alleles genome-wide')

# read multivariate eGWAS results and concatenate all together
resultsFilePath = '/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/'
all_files = glob.glob(os.path.join(resultsFilePath, "envGWAS_results*"))
resultsFile = pd.concat((pd.read_csv(f) for f in all_files), ignore_index = True)
top_hits_coord = pd.read_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/GEA_lead_SNPs_list_with_coord.csv')
sig_hits_coord = pd.read_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/GEA_significant_SNPs_list_with_coord.csv')

## do not filter SNPs for new envGWAS with 7 climate variables
## import maf-filtered SNPS (maf > 0.05) generated from util_functions::read_MAF_file()
# maf_filter = pd.read_table('/group/runciegrp2/Projects/SeeD//data/maf_filtered_SNPs.txt', 
#                            sep = ' ', 
#                            index_col = None, 
#                            names = ['snp', 'pval'],
#                           skiprows = 1)


## define parameters

bp = 25000 
# pval = 1e-2
# look at all SNPs, irrelevant of p-value
p1_SNPs = resultsFile
# p1_SNPs = p1_SNPs[p1_SNPs['snp'].isin(maf_filter['snp'])]

## calculate inv4m
inv4m_coord = top_hits_coord.loc[16]
## calculate hsftf9
hsftf9_coord = top_hits_coord.loc[30]
## drop inv4m from this
# top_hits_coord.drop(16, inplace = True)

# for chrom in range(1, 11):
#     genotypeFile = pd.DataFrame()
#     chr_p1 = p1_SNPs[p1_SNPs['Chr'] == chrom]
#     filtered_chr_p1 = pd.DataFrame()
#     for top_pos in top_hits_coord[top_hits_coord['CHR'] == chrom]['BP']:
#         filtered_chr_p1 = pd.concat([filtered_chr_p1, chr_p1[abs(chr_p1['pos'] - top_pos) < bp]])

#     if filtered_chr_p1.shape[0] == 0:
#         continue
#     chr_genotypeFile = pd.read_csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome/chr%02d.012.csv' % chrom,
#         header = 0,
#         usecols = filtered_chr_p1['snp'],
#         dtype = np.float64,
#         engine = 'c')
#     print('Loaded genotypeFile chromosome %d' % chrom)
#     # genotypeFile = pd.concat([genotypeFile, chr_genotypeFile], axis = 1)
#     # print('Concatenated genotypeFile chromosome %d' % chrom)
#     chr_genotypeFile.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_SNPs_within_significant_windows_chr%02d.012.csv' % chrom, index = False)
#     print('Printed genotypeFile')

#     freqcorr = np.square(chr_genotypeFile.corr())
#     print('Calculated correlation matrix')
#     freqcorr.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_local_corr_matrix_chr%02d.csv' % chrom, index = True)
#     print('Printed correlation matrix')


## use only inv4m 
# for chrom in range(4, 5):
#     genotypeFile = pd.DataFrame()
#     chr_p1 = p1_SNPs[p1_SNPs['Chr'] == chrom]
#     filtered_chr_p1 = chr_p1[abs(chr_p1['pos'] - inv4m_coord['BP']) < 600000]

#     chr_genotypeFile = pd.read_csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome/chr%02d.012.csv' % chrom,
#         header = 0,
#         usecols = filtered_chr_p1['snp'],
#         dtype = np.float64,
#         engine = 'c')
#     print('Loaded genotypeFile chromosome %d' % chrom)
#     # genotypeFile = pd.concat([genotypeFile, chr_genotypeFile], axis = 1)
#     # print('Concatenated genotypeFile chromosome %d' % chrom)
#     chr_genotypeFile.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_SNPs_within_significant_windows_inv4m.012.csv', index = False)
#     print('Printed genotypeFile')

#     freqcorr = np.square(chr_genotypeFile.corr())
#     print('Calculated correlation matrix for inv4m')
#     freqcorr.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_local_corr_matrix_inv4m.csv', index = True)
#     print('Printed correlation matrix for inv4m')


## use only hsftf9 locus 
for chrom in range(9, 10):
    genotypeFile = pd.DataFrame()
    chr_p1 = p1_SNPs[p1_SNPs['Chr'] == chrom]
    # filtered_chr_p1 = chr_p1[abs(chr_p1['pos'] - hsftf9_coord['BP']) < 600000]

    chr_genotypeFile = pd.read_csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome/chr%02d.012.csv' % chrom,
        header = 0,
        usecols = ['S9_108896190', 'S9_110229263', 'S9_110887326', 'S9_111510530', 'S9_148365695'],
        dtype = np.float64,
        engine = 'c')
    print('Loaded genotypeFile chromosome %d' % chrom)
    # genotypeFile = pd.concat([genotypeFile, chr_genotypeFile], axis = 1)
    # print('Concatenated genotypeFile chromosome %d' % chrom)
    chr_genotypeFile.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_SNPs_within_significant_windows_hsftf9.012.csv', index = False)
    print('Printed genotypeFile')

    freqcorr = np.square(chr_genotypeFile.corr())
    print('Calculated correlation matrix for hsftf9')
    freqcorr.to_csv('/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/LD_analysis/LD_local_corr_matrix_hsftf9.csv', index = True)
    print('Printed correlation matrix for hsftf9')


