####################################################################################
# Scripts for candidate gene discovery
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################
source(here::here('config.R'))
{
  v4_coords = read.delim('~/Documents/Projects/gea-adaptation/data/gbs2.7_agpv4.bed', header = F)
  colnames(v4_coords) = c('CHR', 'BP', 'END', 'SNP', 'STRAND')
  v4_coords$CHR = as.integer(v4_coords$CHR)
  v4_coords = v4_coords %>% drop_na(CHR)
  v4_coords$v4_SNP = paste0('S', v4_coords$CHR, '_', v4_coords$END)
}
clumped_SNPs = load_clumped()
clumped_SNPs$clumped = lapply(clumped_SNPs$clumped, function(x) {gsub("\\[|\\]|'|\\s", "", x)})
clumped_SNPs_pivot = clumped_SNPs %>% 
  separate_rows(clumped, sep = ',', convert = TRUE)
# clumped_SNPs_pivot = merge(clumped_SNPs_pivot %>% select(-c(X, num_clumped)), v4_coords[c('SNP', 'CHR', 'END', 'v4_SNP')], by.x = 'clumped', by.y = 'v4_SNP')
colnames(clumped_SNPs_pivot) = c('X', 'lead_SNP', 'Chr', 'BP', 'P', 'num_clumped', 'clumped')
candidate_snps = clumped_SNPs_pivot %>% pull('lead_SNP')
# clumped_snps <<- lapply(genome_clumped_03$clumped, function(x) {str_split(x, ",")})
# clumped_snps = unlist(clumped_snps)
non_lead_snps = clumped_SNPs_pivot %>% pull('clumped')

library(tibble)
v4_candidates = v4_coords %>% filter(v4_SNP %in% candidate_snps | v4_SNP %in% non_lead_snps) %>% select(c('CHR', 'BP', 'END', 'v4_SNP', 'STRAND'))
v4_candidates = add_column(v4_candidates, score = rep(0, nrow(v4_candidates)), .after = 4)
colnames(v4_candidates) = c('chrom', 'start', 'end', 'name', 'score', 'strand')
v4_candidates = v4_candidates[order(v4_candidates$chrom, v4_candidates$start),]
# write.table(v4_candidates, here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'GEA_candidate_genes_v4.bed'), quote = F, sep = '\t', row.names = F)

#egrep '^[1-9]' ~/Documents/Projects/gea-adaptation/data/Zea_mays.AGPv4.36.gff3 | grep "gene" | grep -v "transcript" | sort -k1,1n -k4,4n | bedtools closest -D a -t first -a GEA_candidate_genes_v4.bed -b stdin > GEA_candidate_genes_closest.txt

closest = fread(here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'GEA_candidate_genes_closest.txt'))
colnames(closest) = c('chrom', 'start', 'end', 'SNP', 'score', 'strand', 
                      'gene_chrom', 'source', 'type', 'gene_start', 'gene_end', 'gene_score', 'gene_strand', 'V14',
                      'gene_description', 'dist_from_SNP')
closest_genes = extract(closest, gene_description, into = 'gene_ID', regex = "(?<=gene_id=)(.*)(?=;logic_name)", remove = F)

clumped_genes = merge(clumped_SNPs_pivot, closest_genes, by.x = 'clumped', by.y = 'SNP')
ggplot(clumped_genes %>% 
         group_by(lead_SNP) %>% 
         # filter(!lead_SNP %in% c('S4_181163723', 'S2_229299879', 'S4_172772017', 'S4_170970552')) %>% 
         summarise(across(gene_ID, n_distinct, na.rm = T)),
       aes(x = gene_ID)) +
  geom_histogram(bins = 7) +
  ggtitle("How many clumped peaks hit one gene?") +
  xlab('# unique gene hits in peak') +
  ylab('# clumped peaks')
  

ggplot(clumped_genes %>% 
         aes(x = gene_ID))

closest_genes
View(closest_genes)
