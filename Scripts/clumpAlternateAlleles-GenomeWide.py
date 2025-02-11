import numpy as np
import pandas as pd
import glob
import os

print('Start clumping alleles genome-wide')



# read multivariate eGWAS results and concatenate all together
resultsFilePath = '/group/runciegrp2/Projects/SeeD/Analyses/GEA_output/multivariate_results_unstructured'
all_files = glob.glob(os.path.join(resultsFilePath, "envGWAS_results*"))
resultsFile = pd.concat((pd.read_csv(f) for f in all_files), ignore_index = True)

## do not filter SNPs for new envGWAS with 7 climate variables
## import maf-filtered SNPS (maf > 0.05) generated from util_functions::read_MAF_file()
# maf_filter = pd.read_table('/group/runciegrp2/Projects/SeeD//data/maf_filtered_SNPs.txt', 
#                            sep = ' ', 
#                            index_col = None, 
#                            names = ['snp', 'pval'],
#                           skiprows = 1)


## define parameters
# chrom = 4
r2 = 0.3
bp = 15000000 # only used for window-based clumping, not used in this script
pval = 1e-250
# only look at SNPs below certain p-value threshold
p1_SNPs = resultsFile[resultsFile['X..Pvalue'] < pval]
# p1_SNPs = p1_SNPs[p1_SNPs['snp'].isin(maf_filter['snp'])]


genotypeFile = pd.DataFrame()
for chrom in range(1, 11):
    
    chr_p1 = p1_SNPs[p1_SNPs['Chr'] == chrom]
    chr_genotypeFile = pd.read_csv('/group/runciegrp2/Projects/SeeD/Genetic_data/Imputed_V4/genotypes_by_chromosome/chr%02d.012.csv' % chrom,
        header = 0,
        usecols = chr_p1['snp'],
        dtype = np.float64,
        engine = 'c')
    print('Loaded genotypeFile chromosome %d' % chrom)
    genotypeFile = pd.concat([genotypeFile, chr_genotypeFile], axis = 1)
    print('Concatenated genotypeFile chromosome %d' % chrom)

genotypeFile.to_csv('%s/clumped/significant-SNPs.012.csv' % resultsFilePath, index = False)
print('Printed genotypeFile', flush = True)

# genotypeFile = pd.read_csv('%s/clumped/significant-SNPs.012.csv' % resultsFilePath)


freqcorr = np.square(genotypeFile.corr())
print('Calculated correlation matrix', flush = True)
freqcorr.to_csv('%s/clumped/corr_matrix.csv' % resultsFilePath, index = True)
print('Printed correlation matrix', flush = True)

# freqcorr = pd.read_csv('%s/clumped/corr_matrix.csv' % resultsFilePath)

freqcorrmelt = pd.melt(freqcorr.reset_index(), 
        id_vars = 'index', 
        value_vars = freqcorr.columns,
        var_name = 'B',
        value_name = 'R2')
freqcorrmelt['R2'] = pd.to_numeric(freqcorrmelt['R2'])
print('Melted correlation matrix', flush = True)



# sort by p-value, greedy algorithm gives priority to most significant
p1_SNPs.sort_values(by = 'X..Pvalue', inplace = True)
p1_df = pd.DataFrame(columns = ['SNP', 'Chr', 'BP', 'P', 'num_clumped', 'clumped'])
test = p1_SNPs

while not p1_SNPs.empty:
    row = p1_SNPs.iloc[0]
    snp = row['snp']
    pos = row['pos']
    this_pval = row['X..Pvalue']
    chrom = row['Chr']
    print(snp)
    # calculate pairwise R^2 value between target SNP and nearby SNPs
    clumped_snps = freqcorrmelt[(freqcorrmelt['index'] == snp) & (freqcorrmelt['R2'] > r2)]
    num_clumped_snps = len(clumped_snps['index'])
    p1_df = pd.concat([p1_df, 
        pd.DataFrame({'Chr': [chrom], 
                       'SNP': [snp], 
                       'BP': [pos], 
                       'P': [this_pval], 
                        'num_clumped': [num_clumped_snps],
                       'clumped': [str(clumped_snps['B'].to_list())]
                       })],
                     ignore_index = True)
    # pop goes the weasel
    p1_SNPs.drop(p1_SNPs[p1_SNPs['snp'].isin(list(clumped_snps['B']))].index, inplace = True)
    # p1_SNPs = p1_SNPs.iloc[1:,:]

p1_df.to_csv('%s/clumped/envGWAS_results.genomeclumped_%s.csv' % (resultsFilePath, str(pval)))
print('Finished clumping', flush = True)



