# Bar plots of GSEA analysis results

# CREATED : 3/18/2020 by Cherry Gao
#           Adapted from bar_plot_gene_sets.R
# UPDATED : 3/30/2020 
#           Order the bars in a logical way.


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
gsea_results_dir = paste0('gsea/gsea_results/fdr_',deseq_alpha) # where GSEA results should be called from
plot_save_dir = paste0('gsea/gsea_results/fdr_',deseq_alpha,'/bar_plots_sets/ordered_by_FDR/path_label_bars') # where to save the bar plots
#   dir.create(plot_save_dir)

# define the pairwise comparisons range
pairwise_range = c(2,3) # c(2:length(cond_pair_list))

gsea_kegg_for_gene_bar = list() # initialize container of KEGG to make bar plots of
# load each pairwise GSEA result and make bar graph
for (z in pairwise_range){ # skip 0m vs. 0c
  
  pairwise_title = paste0(cond_pair_list[[z]][1],'_vs_',cond_pair_list[[z]][2])
  fn = paste0(pairwise_title,'_gsea_results_clean.Rdata')
  load(file=paste0(gsea_results_dir,'/',fn))
  gsea_kegg_for_gene_bar[[z]] = matrix()

  # reassign loaded variables for simplicity
  neg = neg_genesets_table 
  pos = pos_genesets_table

  # add "negative" or "positive" factor to data according to NES
  neg['enrichment'] = 'negative'
  pos['enrichment'] = 'positive'
  
  # FDR threshold: extract indexes in data_concat of significant gene sets to be plotted as bars 
  i_pos =  which(pos$FDR.q.val < gsea_fdr_thresh) # extract indexes in data_concat of genes to be plotted as bars
  i_neg =  which(neg$FDR.q.val < gsea_fdr_thresh) # extract indexes in data_concat of genes to be plotted as bars
  
  # subset on those that pass FDR threshold
  pos_sub = pos[i_pos,]
  neg_sub = neg[i_neg,]
  
  # order by FDR q-value
  pos_sub = pos_sub[order(pos_sub$FDR.q.val),]
  neg_sub = neg_sub[order(neg_sub$FDR.q.val),]
  
  # prepare data for plot, order positive and negative values separately by significance
  data_combined = rbind(pos_sub,neg_sub)

  # order of the bars should be order of the indexes
  data_combined['order'] = rev(c(1:length(data_combined[,1])))
 #  data_combined['path_labe'] = unlist(strsplit(as.character(data_combined$NAME),'_'))
  
  # string split and convert to lower case for bar labels
  split_one = str_split(data_combined$NAME,'_',simplify=T)[,2]
  split_path_only = str_split(split_one," \\[PATH:VCT",simplify=TRUE)[,1]
  data_combined['bar_label'] = paste0(tolower(split_path_only),' (',data_combined$SIZE,')')
  
  # ggplot : bar graph of all gene sets
  title_text = paste(pairwise_title,'\n GSEA results; KEGG PATHWAYs as genesets', '\n FDR threshold = ',gsea_fdr_thresh, '\n ordered by FDR q-values (black labels), small to large')
  
  p = ggplot(data_combined, aes(x=order, y=NES, fill=enrichment)) + geom_col(width = 0.9) +
      theme_bw() + # get rid of grey background
      scale_fill_manual(values = c('#969696', '#6baed6') ) + # false, true bar colors
   #   ylim(-2.5,3) + # set x limit the same for all pairwise comparisons (ylim because plot is flipped)
      ylim(-2.5,4) + # for 10m vs 10c and 60m vs 60c
      coord_flip() + 
      scale_x_continuous(breaks=NULL) # get rid of vertical tick lines
  
  # label bars with pathway name on the y axis
#  p = p + scale_x_continuous(
#            breaks = data_combined$order,
#            labels = data_combined$NAME,
#            expand = c(0,0)
#            ) + coord_flip()
  
  # plot labels
  p = p + ggtitle(title_text) + ylab('Normalized Enrichment Score') + xlab('KEGG pathway gene sets') +
      #    geom_text(aes(label=FDR.q.val, y=0.5*NES), colour="black", size=1) # FDR q-value label on bars
         geom_text(aes(label=bar_label, y=0.5*NES), colour="black", size=1) # pathway names on bars
    
  p = p + theme(text = element_text(size=4),
                plot.title = element_text(hjust = 1, size = 4, color='black')) 
  
  p # show plot if not in for-loop
  
  # save figure (eps ends up in broken text; pdf is better)
  save_fn = paste0(plot_save_dir,'/',pairwise_title,'_gsea_results_bar_graph_v2.pdf')
  ggsave(save_fn, dpi=300, dev='pdf', height = 3, width = 4, units = "in") 
  # ggsave(save_fn, dpi=300, dev='png', height = 3, width = 4, units = "in") 
  
  # extract KEGG pathway codes
  for (path in 1:length(data_combined$NAME)){
    path_name = as.character(data_combined$NAME[path])
    gsea_kegg_for_gene_bar[[z]] = append(gsea_kegg_for_gene_bar[[z]], unlist(strsplit(path_name, '_'))[1])
  }
}

# save the list of KEGG pathway codes to plot next
# save(gsea_kegg_for_gene_bar, file=paste0(gsea_results_dir,'/GSEA_results_kegg_path_genes_to_plot.Rdata'))
  
  