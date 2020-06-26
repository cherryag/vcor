# Bar plots of GSEA analysis results

# CREATED : 3/18/2020 by Cherry Gao
#           Adapted from bar_plot_gene_sets.R


#---------SUMMARY---------
# 1) Load each pairwise GSEA results.
# 2) Cutoff by FDR q value threshold.
# 3) Make and save bar graphs of KEGG pathways that are positively and negatively enriched.
# 4) Extract the KEGG pathway codes for next step: bar plot of genes in each KEGG pathway.
#-------------------------

# libraries
library(ggplot2)
library(stringr)
library(tidyverse)


############### SET UP ###############
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load('pairwise_list.Rdata')
deseq_alpha = 0.01
gsea_fdr_thresh = 0.25 # 25% was recommended by GSEA

# GSEA results folder
gsea_results_dir = paste0('gsea/gsea_results/fdr_',deseq_alpha)

gsea_kegg_for_gene_bar = list() # initialize container of KEGG to make bar plots of
# load each pairwise GSEA result and make bar graph
for (z in 2:length(cond_pair_list)){ # skip 0m vs. 0c
  
  pairwise_title = paste0(cond_pair_list[[z]][1],'_vs_',cond_pair_list[[z]][2])
  fn = paste0(pairwise_title,'_gsea_results_clean.Rdata')
  load(file=paste0(gsea_results_dir,'/',fn))
  gsea_kegg_for_gene_bar[[z]] = matrix()

  # reassign loaded variables for simplicity
  neg = neg_genesets_table 
  pos = pos_genesets_table

  # FDR threshold: extract indexes in data_concat of significant gene sets to be plotted as bars 
  i_pos =  which(pos$FDR.q.val < gsea_fdr_thresh) # extract indexes in data_concat of genes to be plotted as bars
  i_neg =  which(neg$FDR.q.val < gsea_fdr_thresh) # extract indexes in data_concat of genes to be plotted as bars
  
  # prepare data for plot
  data_combined = rbind(pos[i_pos,], neg[i_neg,])
   
  # ggplot : bar graph of all gene sets
  title_text = paste(pairwise_title,'\n GSEA results; KEGG PATHWAYs as genesets', '\n FDR threshold = ',gsea_fdr_thresh)
  
  p = ggplot(data_combined, aes(x = NAME, y = NES, fill=(NES > 0))) + 
        theme_bw() + # get rid of grey background
        geom_col(width = 0.9) + # bar width with respect to each other (1 = bars touch each other)
        coord_flip() + # horizontal bars
        scale_fill_manual(values = c('purple', 'orange') )  # false, true bar colors
    
  # labels
  p = p + ggtitle(title_text) + ylab('Normalized Enrichment Score') + xlab('KEGG pathway gene sets')
  p = p + theme(text = element_text(size=4),
                plot.title = element_text(hjust = 1, size = 4, color='black'))
  
  # p # if not in for-loop
  #print(p)
  
  # save figure
  save_fn = paste0(gsea_results_dir,'/',pairwise_title,'_gsea_results_bar_graph.eps')
  ggsave(save_fn, dpi=300, dev='eps', height = 3, width = 4, units = "in") 
  ggsave(save_fn, dpi=300, dev='png', height = 3, width = 4, units = "in") 
  
  # extract KEGG pathway codes
  for (path in 1:length(data_combined$NAME)){
    path_name = as.character(data_combined$NAME[path])
    gsea_kegg_for_gene_bar[[z]] = append(gsea_kegg_for_gene_bar[[z]], unlist(strsplit(path_name, '_'))[1])
  }
}

# save the list of KEGG pathway codes to plot next
save(gsea_kegg_for_gene_bar, file=paste0(gsea_results_dir,'/GSEA_results_kegg_path_genes_to_plot.Rdata'))
  
  