# Bar plots of log2FC of sets of genes, color-coded by padj
# Adapted from old script named log2plot.R from 2016


#---------- HISTORY ----------
# CREATED : 1/30/2020 by Cherry Gao
# UPDATED : 2/9/2020  added virulence gene set (manual)
# UPDATED : 3/18/2020
#           Completely re-vamped, now referencing JGI gene ID#s in KEGG pathways
#-----------------------------

#----------SUMMARY----------
# 1) Identify the KEGG Pathway whose individual genes are to be plotted
# 2) Identify the JGI gene ID and indexes of those genes in data_concat (DESeq2 result file)
# 3) Create and save bar plot for each set of genes. 
# 4) Cycle through each pairwise comparison. 
#---------------------------

# libraries
library(ggplot2)
library(stringr)



#---------- SET UP
rm(list = ls()) # clear environment
deseq_alpha = 0.01
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load('pairwise_list.Rdata')
load('dataset/kegg_downloaded/kegg_pathways_jgi.Rdata')
load(file=paste0('gsea/gsea_results/fdr_',deseq_alpha,'/GSEA_results_kegg_path_genes_to_plot.Rdata')) # list of KEGG pathways from 11_bar_plot_gsea.R


# GSEA results folder
gsea_results_dir = paste0('gsea/gsea_results/fdr_',deseq_alpha)
deseq_results_dir = 'deseq/deseq_kegg_concat'

# padj (DESeq2 result) threshold for coloring bars
padj_thresh = 0.01 


for (z in 2:length(cond_pair_list)){ # FOR EACH PAIRWISE GSEA, LOAD data_concat (skip 0m vs. 0c)
  
  pairwise_title = paste0(cond_pair_list[[z]][1],'_vs_',cond_pair_list[[z]][2])
  fn = paste0(pairwise_title,'_alpha_',deseq_alpha,'_KEGGconcat.Rdata') # load data_concat
  load(file=paste0(deseq_results_dir,'/',fn))
  
  # housekeeping on data_concat
  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
  data_concat$padj = as.numeric(paste(data_concat$padj))

  # list of KEGG pathways to look at for this pairwise comparison
  kegg_pathway_num = gsea_kegg_for_gene_bar[[z]]
  kegg_pathway_num = kegg_pathway_num[!is.na(kegg_pathway_num)] # remove the first NA

  for (j in 1:length(kegg_pathway_num)){ # FOR EACH KEGG PATHWAY, make a bar plot of all its genes and log2FC

    data_subset = list() # initialize
      
    index_kegg_path = which(names(kegg_list_jgi) == kegg_pathway_num[j]) # extract index of the pathway in KEGG Pathway database

    # sequential processing of JGI gene IDs in the KEGG pathway
    path_jgi_gene_list = kegg_list_jgi[index_kegg_path]
    path_jgi_gene_list = path_jgi_gene_list[[1]]
      path_description = path_jgi_gene_list[1] # extrat pathway description
    path_jgi_gene_list = as.numeric(path_jgi_gene_list[2:length(path_jgi_gene_list)]) # only keep the JGI gene IDs; convert to integer for matching

    # find all indexes in data_concat of KEGG pathway genes
    deseq_jgi_gene_list = data_concat$seq_data_gene_id # extract all JGI gene IDs from the data_concat
    i_jgi_kegg = match(path_jgi_gene_list, deseq_jgi_gene_list) # extract indexes of JGI gene IDs of interest from data_concat

    # print warning if data_concat does not contain 5022 genes. 
    if (!length(deseq_jgi_gene_list) == 5022) {print(paste('DESeq result dataset does not have 5022 genes (',pairwise_title,')'))} # there should be all genes in each data_concat
    
    # print warning if # of JGI genes between KEGG pathway database and DESeq results do not match
    if (!length(path_jgi_gene_list) == length(i_jgi_kegg)) {print(paste('mismatch in number of JGI genes between KEGG pathway database and those found in data_concat (',pairwise_title,')'))} 

    # extract DESeq result data of the genes to plot
    data_subset = data_concat[i_jgi_kegg,]
    data_subset$bar_labels = paste0(data_concat$geneProduct[i_jgi_kegg],"_", data_concat$jgi_gene_id[i_jgi_kegg]) # for bar labels - concatenate JGI gene ID#

    # ------- create bar plot
    title_text = paste(pairwise_title, '\n', path_description,
                       '\n Log2 fold change and significance of genes in KEGG pathways from GSEA results', 
                       '\n padj threshold = ',padj_thresh)
    
    p <- ggplot(data_subset, aes(x = foldchange, y = bar_labels, fill=(padj < padj_thresh)), stat='identity') +
          geom_col(width = 0.9)  + 
          theme_bw() +
          # coord_flip() +
          scale_fill_manual(values = c('grey', 'red') ) #+  # false, true bar colors
       #   xlim(-5, 5.5)

    # labels
    p = p + ggtitle(title_text) + ylab('Log2 Fold Change') + xlab('gene product _ JGI gene ID')
    p = p + theme(text = element_text(size=4),
                  plot.title = element_text(hjust = 1, size = 4, color='black'))
    
   # p # if not in for-loop
   # print(p) # if in for-loop

    # save figure (png first, for inspection)
    save_fn = paste0(gsea_results_dir,'/bar_plots_genes/',pairwise_title,'_kegg_genes_bar_',kegg_pathway_num[j])
#    ggsave(paste0(save_fn,'.eps'), dpi=300, dev='eps', height = 3, width = 4, units = "in") 
    ggsave(paste0(save_fn,'.png'), dpi=300, dev='png', height = 6, width = 6, units = "in") 

    } # for each KEGG pathway
} # for each pairwise comparison






  