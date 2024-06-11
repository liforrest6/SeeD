####################################################################################
# Scripts for candidate gene discovery
#
# Author: Forrest Li
# Functions for running analyses
####################################################################################
source(here::here('config.R'))
library(gdata)
{
  v4_coords = read.delim('~/Documents/Projects/SeeD/Genetic_data/gbs2.7_agpv4.bed', header = F)
  colnames(v4_coords) = c('CHR', 'BP', 'END', 'SNP', 'STRAND')
  v4_coords$CHR = as.integer(v4_coords$CHR)
  v4_coords = v4_coords %>% drop_na(CHR)
  v4_coords$v4_SNP = paste0('S', v4_coords$CHR, '_', v4_coords$END)
}
clumped_SNPs = load_clumped()
top_hits_clumped = load_clumped()$SNP
clumped_SNPs$clumped = lapply(clumped_SNPs$clumped, function(x) {gsub("\\[|\\]|'|\\s", "", x)})
clumped_SNPs_pivot = clumped_SNPs %>% 
  separate_rows(clumped, sep = ',', convert = TRUE)
# clumped_SNPs_pivot = merge(clumped_SNPs_pivot %>% select(-c(X, num_clumped)), v4_coords[c('SNP', 'CHR', 'END', 'v4_SNP')], by.x = 'clumped', by.y = 'v4_SNP')
colnames(clumped_SNPs_pivot) = c('X', 'lead_SNP', 'Chr', 'BP', 'P', 'num_clumped', 'clumped')
candidate_snps = clumped_SNPs_pivot %>% pull('lead_SNP')
# clumped_snps <<- lapply(genome_clumped_03$clumped, function(x) {str_split(x, ",")})
# clumped_snps = unlist(clumped_snps)
non_lead_snps = clumped_SNPs_pivot %>% pull('clumped')

v4_candidates = v4_coords %>% filter(v4_SNP %in% candidate_snps | v4_SNP %in% non_lead_snps) %>% dplyr::select(c('CHR', 'BP', 'END', 'v4_SNP', 'STRAND'))
v4_candidates = add_column(v4_candidates, score = rep(0, nrow(v4_candidates)), .after = 4)
colnames(v4_candidates) = c('chrom', 'start', 'end', 'name', 'score', 'strand')
v4_candidates = v4_candidates[order(v4_candidates$chrom, v4_candidates$start),]
write.table(v4_candidates, here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'GEA_candidate_SNP_coordinates_v4.bed'), quote = F, sep = '\t', row.names = F, col.names = F)

## make sure to re-sort bed file to lexicographic because I hate bedtools
##  gunzip -c ~/Documents/Data/Zea_mays.AGPv4.36.gff3 | egrep '^[1-9]' - | grep "gene" | grep -v "transcript" | sort -k1,1n -k4,4n | \
## bedtools closest -D a -t first -a GEA_candidate_SNP_coordinates_v4.bed -b stdin > ~/Documents/Projects/SeeD/Analyses/GEA_output/multivariate_results/clumped/GEA_candidate_genes_closest.txt

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
  geom_histogram(bins = 100) +
  ggtitle("How many clumped peaks hit one gene?") +
  xlab('# unique gene hits in peak') +
  ylab('# clumped peaks')

genes_all = read.table('~/Documents/Data/zeamays_genemodels/genes_all.txt', sep = '\t', fill = T, header = T)
genes_classical = read.table('~/Documents/Data/zeamays_genemodels/genes_classical.txt', sep = '\t', fill = T, header = T, quote = '')
genes_maizegdb = read.table('~/Documents/Data/zeamays_genemodels/genes_maizegdb.txt', sep = '\t', fill = T, header = T, quote = '')

all_sourced_models = genes_all %>% filter(v4.Gene.Model.ID %in% clumped_genes$gene_ID)
## we see 89 genes with any kind of model
genes_classical %>% filter(v4.Gene.Model.ID %in% closest_genes$gene_ID)
## we don't see any classical genes that have any of the V4 IDs
genes_maizegdb %>% filter(v4.Gene.Model.ID %in% closest_genes$gene_ID)
## we see 16 maizegdb genes
createLink = function(gene) {
  return(sprintf('https://www.maizegdb.org/gene_center/gene/%s', gene))
}

gene_links = lapply(unique(all_sourced_models %>% filter(Full.Name != '') %>% pull(v4.Gene.Model.ID)), createLink)

## manually selected from all_sourced_models
final_genes = c('Zm00001d052227', 'Zm00001d052060', 'Zm00001d052063',
                'Zm00001d052185', 'Zm00001d052197', 'Zm00001d051847',
                'Zm00001d051902', 'Zm00001d051945', 'Zm00001d051976',
                'Zm00001d052018', 'Zm00001d052269', 'Zm00001d048041')

final_models = all_sourced_models %>% filter(v4.Gene.Model.ID %in% final_genes)
final_coordinates = merge(final_models[c('v4.Gene.Model.ID', 'Gene.Symbol', 'Full.Name')], closest_genes[c('gene_ID', 'gene_chrom', 'gene_start', 'gene_end')], by.x = 'v4.Gene.Model.ID', by.y = 'gene_ID') %>% distinct()
## 10 lead SNPs that have clumped SNPs that have a closest gene that has a sourced gene model
clumped_genes %>% filter(gene_ID %in% all_sourced_models$v4.Gene.Model.ID) %>% pull(lead_SNP) %>% unique()
## but only Inv4m and hsftf9 have lead SNPs from manually selected candidate genes
clumped_genes %>% filter(gene_ID %in% final_models$v4.Gene.Model.ID) %>% pull(lead_SNP) %>% unique()

## look for candidate genes that are closest to only lead SNPs
lead_candidate_genes = clumped_genes %>% filter(clumped %in% candidate_snps) %>% pull(gene_ID)
## there are only 9 sourced models corresponding to those lead SNPs
lead_models = all_sourced_models %>% filter(v4.Gene.Model.ID %in% lead_candidate_genes)
lead_models_clumped_genes = clumped_genes %>% filter(gene_ID %in% lead_models$v4.Gene.Model.ID) %>% dplyr::select(clumped, lead_SNP, gene_chrom, gene_start, gene_ID)
clumped_genes %>% filter(lead_SNP %in% lead_models_clumped_genes$lead_SNP)
clumped_genes %>% filter(lead_SNP %in% candidate_snps)

gene_LD_plots = list()
nearby_SNPs = c()
# SNP_start_in_candidate_genes = closest_genes %>% filter(gene_ID %in% final_genes) %>% pull(gene_start)
# SNP_end_in_candidate_genes = closest_genes %>% filter(gene_ID %in% final_genes) %>% pull(gene_end)
for(i in 1:12){
  SNP_start_in_candidate_genes = final_coordinates[i, 'gene_start']
  SNP_end_in_candidate_genes = final_coordinates[i, 'gene_end']
  nearby_SNPs = gea_results %>% filter(CHR == final_coordinates[i, 'gene_chrom'] & BP < SNP_end_in_candidate_genes + 5000 & BP > SNP_start_in_candidate_genes - 5000) %>% 
    mutate(BP_centered = BP - SNP_start_in_candidate_genes)
  plot_i = ggplot(nearby_SNPs, aes(x = BP_centered, y = -log(P))) +
    geom_point() +
    geom_vline(xintercept = 0, color = 'red') +
    geom_vline(xintercept = SNP_end_in_candidate_genes - SNP_start_in_candidate_genes, color = 'red') +
    # xlim(c(SNP_start_in_candidate_genes - 5000, SNP_end_in_candidate_genes + 5000)) +
    xlim(c(-5000, SNP_end_in_candidate_genes - SNP_start_in_candidate_genes + 5000)) +
    ylim(0, 20) +
    ggtitle(paste(final_coordinates[i, 'Gene.Symbol'], final_coordinates[i, 'v4.Gene.Model.ID'], sep = ', '))
  
  gene_LD_plots[[i]] = plot_i
}

for(i in 1:12){
  SNP_start_in_candidate_genes = final_coordinates[i, 'gene_start']
  SNP_end_in_candidate_genes = final_coordinates[i, 'gene_end']
  nearby_SNPs = gea_results %>% filter(CHR == final_coordinates[i, 'gene_chrom'] & BP < SNP_end_in_candidate_genes + 5000 & BP > SNP_start_in_candidate_genes - 5000)
  plot_i = ggplot(nearby_SNPs, aes(x = BP, y = -log(P))) +
    geom_point() +
    geom_vline(xintercept = SNP_start_in_candidate_genes, color = 'red') +
    geom_vline(xintercept = SNP_end_in_candidate_genes, color = 'red') +
    xlim(c(SNP_start_in_candidate_genes - 5000, SNP_end_in_candidate_genes + 5000)) +
    ylim(0, 20) +
    ggtitle(paste(final_coordinates[i, 'Gene.Symbol'], final_coordinates[i, 'v4.Gene.Model.ID'], sep = ', '))
  
  gene_LD_plots[[i]] = plot_i
}


candidate_gene_LD = plot_grid(plotlist = gene_LD_plots, nrow = 3, ncol = 4)
png(here(plot_dir, 'Manuscript', 'candidate_gene_LD.png'), width = 491, height = 491)
print(candidate_gene_LD)
dev.off()

closest_genes
View(closest_genes)

### LD analysis - decay plot based on r^2 scores of clumped significant genes

corr_matrix = read.csv(here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'corr_matrix.csv'))
corr_matrix_melt = pivot_longer(corr_matrix, cols = starts_with('S'), names_to = 'Clumped SNP', values_to = 'r2')
corr_matrix_melt = merge(corr_matrix_melt, v4_coords[c('v4_SNP', 'BP', 'CHR')], by.x = 'Clumped SNP', by.y = 'v4_SNP')

LD_plots = lapply(top_hits_clumped, function(lead_SNP) {
  lead_SNP_pos = v4_coords[v4_coords$v4_SNP == lead_SNP, 'BP']
  lead_SNP_chr = v4_coords[v4_coords$v4_SNP == lead_SNP, 'CHR']
  lead_corrs = corr_matrix_melt %>% filter(X == lead_SNP, CHR == lead_SNP_chr)
  
  lead_corrs = lead_corrs %>% filter(r2 > 0.5)
  if(nrow(lead_corrs) > 1){
    ggplot(lead_corrs, aes(x = BP - lead_SNP_pos, y = r2)) + 
      geom_point() +
      geom_vline(xintercept = 0, color = 'red') +
      ggtitle(lead_SNP) +
      ylim(c(0.5, 1))
  }

})




# write.csv(v4_coords %>% filter(v4_SNP %in% candidate_snps), 
#           here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'GEA_lead_SNPs_list_with_coord.csv'),
#           row.names = F, quote = F)
# write.csv(v4_coords %>% filter(v4_SNP %in% non_lead_snps), 
#           here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'GEA_significant_SNPs_list_with_coord.csv'),
#           row.names = F, quote = F)


compile_correlations = function(chr){
  corr_matrix = data.frame(fread(here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'LD_analysis', sprintf('LD_local_corr_matrix_chr%02d.csv', chr))))
  corr_matrix_melt = pivot_longer(corr_matrix, cols = starts_with('S'), names_to = 'SNP B', values_to = 'r2')
  corr_matrix_melt
}

all_corr_matrix = data.table::rbindlist(lapply(c(1, 2, 3, 4, 5, 7, 8, 9, 10), compile_correlations))
all_corr_matrix = merge(all_corr_matrix, v4_coords[c('v4_SNP', 'BP', 'CHR')], by.x = 'SNP B', by.y = 'v4_SNP')
all_corr_matrix = merge(all_corr_matrix, gea_results[c('SNP', 'P')], by.x = 'SNP B', by.y = 'SNP')

r2_threshold = 0
LD_plots = lapply(top_hits_clumped, function(lead_SNP) {
  lead_SNP_pos = v4_coords[v4_coords$v4_SNP == lead_SNP, 'BP']
  lead_SNP_chr = v4_coords[v4_coords$v4_SNP == lead_SNP, 'CHR']
  lead_corrs = all_corr_matrix %>% filter(V1 == lead_SNP, CHR == lead_SNP_chr, 
                                          BP < lead_SNP_pos + 300000, BP > lead_SNP_pos - 300000)
  
  lead_corrs = lead_corrs
  if(nrow(lead_corrs) > 1){
    ggplot(lead_corrs, aes(x = abs(BP - lead_SNP_pos), y = r2, color = P < 1e-5)) + 
      geom_point() +
      geom_vline(xintercept = 0, color = 'green') +
      geom_hline(yintercept = 0.5, color = 'green')+
      scale_color_manual(values = c('black', 'red')) +
      ggtitle(lead_SNP) +
      ylim(c(r2_threshold, 1)) +
      xlim(c(0, 300000)) +
      scale_x_continuous(limits = c(0, 300000), labels = function(x){x / 1000}) +
      xlab('')+
      ylab('')+
      theme(legend.position = 'none')
  }
  
})

## plotting for inv4m alone
{
  inv4m_matrix = data.frame(fread(here(analyses_dir, 'GEA_output', 'multivariate_results', 'clumped', 'LD_analysis', 'LD_local_corr_matrix_inv4m.csv')))
  inv4m_matrix_melt = pivot_longer(inv4m_matrix, cols = starts_with('S'), names_to = 'SNP B', values_to = 'r2')
  inv4m_corr_matrix = merge(inv4m_matrix_melt, v4_coords[c('v4_SNP', 'BP', 'CHR')], by.x = 'SNP B', by.y = 'v4_SNP')
  inv4m_corr_matrix = merge(inv4m_corr_matrix, gea_results[c('SNP', 'P')], by.x = 'SNP B', by.y = 'SNP')
  inv4m_lead_pos = v4_coords[v4_coords$v4_SNP == 'S4_177835031', 'BP']
  
  inv4m_plot = ggplot(inv4m_corr_matrix %>% filter(V1 == 'S4_177835031'), aes(x = abs(BP - inv4m_lead_pos), y = r2, color = P < 1e-5)) + 
    geom_point() +
    geom_vline(xintercept = 0, color = 'green') +
    geom_hline(yintercept = 0.5, color = 'green')+
    scale_color_manual(values = c('black', 'red')) +
    ggtitle('S4_177835031') +
    ylim(c(0, 1)) +
    xlim(c(1, 600000)) +
    scale_x_continuous(limits = c(0, 600000), labels = function(x){x / 1000}) +
    xlab('')+
    ylab('')+
    theme(legend.position = 'none')
}

LD_plots[[2]] = inv4m_plot
all_plots = plot_grid(plotlist = LD_plots[!sapply(LD_plots, is.null)], ncol = 4)


png(here(plot_dir, 'Manuscript', 'candidate_gene_LD.png'), width = 750, height = 1000)
annotate_figure(all_plots, left = textGrob("LD (r2)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Distance from lead SNP (kb)", gp = gpar(cex = 1.3)))
dev.off()

## check if hsftf9 is in LD with inv9f

all_corr_matrix %>% filter(`SNP B` == 'S9_148365695' & )

lead_hsftf9 = 'S9_148365695'
lead_hsftf9_pos = v4_coords[v4_coords$v4_SNP == lead_hsftf9, 'BP']
lead_hsftf9_chr = v4_coords[v4_coords$v4_SNP == lead_hsftf9, 'CHR']
lead_corrs = all_corr_matrix %>% filter(V1 == lead_hsftf9, CHR == lead_hsftf9_chr, 
                                        BP < lead_hsftf9_pos + 300000, BP > lead_hsftf9_pos - 300000)



