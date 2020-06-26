# Bar plots of log2FC of sets of genes, color-coded by padj
# Version: take in list of JGI gene lists for bar graph (particularly, Kimes et al. virulence factors)


#---------- HISTORY ----------
# CREATED : 1/30/2020 by Cherry Gao
#           Adapted from old script named log2plot.R from 2016
# UPDATED : 2/9/2020  added virulence gene set (manual)
# UPDATED : 3/18/2020
#           Completely re-vamped, now referencing JGI gene ID#s in KEGG pathways
# UPDATED : 3/23/2020 (new version)
#           Gene sets defined by lists of JGI gene ID#s from Kimes et al. potential virulence factors.
# UPDATED : 3/25/2020
#           Now easier to create bar plots with any list of JGI gene ID numbers.
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
deseq_results_dir = 'deseq/deseq_kegg_concat' # to load pairwise result files from 
results_dir = 'virulence_genes' # where to save bar plots
load('pairwise_list.Rdata')
load('dataset/virulence_genes_curated/jgi_gene_id_sets_kimes_qs.Rdata') # load Kimes gene set JGI gene ID list

# DEFINE gene sets [ADDED 3/25/2020]
# Kimes gene sets: 
jgi_gene_id_set = jgi_gene_id_set_kimes # if want to run Kimes et al. gene sets
geneset_list = kimes_geneset_list
title_addendum = '(Kimes et al. gene sets)' # added to plot titles
geneset_range = c(15,17) # c(1:length(jgi_gene_id_set)) # index of gene sets (in geneset_list) to bar plot

# Other gene sets, manually defined by JGI gene ID lists
# jgi_gene_id_set = list()
# jgi_gene_id_set[['Type VI secretion system']] = c(647172191, 647172188, 647172194, 647172193, 647172202, 647172204, 647172182, 647172201)
# geneset_list = names(jgi_gene_id_set)
# title_addendum = 'Guillemette et al. genes identified by homology with V. cholerae T6SS'

# padj (DESeq2 result) threshold for coloring bars
padj_thresh = 0.01 

# list of pairwise comparisons to plot
cond_pair_list_set = c(2:3) # c(2:length(cond_pair_list)) 

for (z in cond_pair_list_set){ # FOR EACH PAIRWISE comparison, LOAD data_concat (skip 0m vs. 0c)
  
  # load data_concat
  pairwise_title = paste0(cond_pair_list[[z]][1],'_vs_',cond_pair_list[[z]][2])
  fn = paste0(pairwise_title,'_alpha_',deseq_alpha,'_KEGGconcat.Rdata') # load data_concat
  load(file=paste0(deseq_results_dir,'/',fn))
  
  # housekeeping on data_concat
  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
  data_concat$padj = as.numeric(paste(data_concat$padj))
  
  # extract all JGI gene IDs from the data_concat (5022 genes)
  deseq_jgi_gene_list = data_concat$seq_data_gene_id 
  
  # run through each Kimes gene sets
  for (j in geneset_range){
    
    # list of JGI gene IDs to plot as bars
    gene_set_description = geneset_list[j] # extract gene set description
    jgi_gene_id = as.numeric(jgi_gene_id_set[[j]])
      
    # find all indexes in data_concat of KEGG pathway genes
    i_jgi_kegg = match(jgi_gene_id, deseq_jgi_gene_list) # extract indexes of JGI gene IDs of interest from data_concat
  
    # print warning if data_concat does not contain 5022 genes. 
    if (!length(deseq_jgi_gene_list) == 5022) {print(paste('DESeq result dataset does not have 5022 genes (',pairwise_title,')'))} # there should be all genes in each data_concat
      
    # print warning if # of JGI genes between list of JGI gene IDs and DESeq results do not match
    if (!length(jgi_gene_id) == length(i_jgi_kegg)) {print(paste('mismatch in number of JGI genes between list of JGI gene IDs and those found in data_concat (',pairwise_title,')'))} 
  
    # extract DESeq result data of the genes to plot
    data_subset = list() # initialize container of log2FC data of genes of interest
    data_subset = data_concat[i_jgi_kegg,]
    data_subset$bar_labels = paste0(data_concat$Locus_tag[i_jgi_kegg],"_",data_concat$geneProduct[i_jgi_kegg],"_", data_concat$jgi_gene_id[i_jgi_kegg]) # for bar labels - concatenate JGI gene ID#
    
      # ------- create bar plot
      title_text = paste(pairwise_title, '\n', paste(gene_set_description, title_addendum),
                         '\n Log2 fold change and significance of genes in gene sets from DESeq results', 
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
      save_fn = paste0(results_dir,'/bar_plots/',pairwise_title,'_Kimes_QS_bar_',gene_set_description)
      # ggsave(paste0(save_fn,'.eps'), dpi=300, dev='eps', height = 3, width = 4, units = "in") 
      ggsave(paste0(save_fn,'.png'), dpi=300, dev='png', height = 6, width = 6, units = "in") 
  
    } # end of Kimes gene sets
} # for each pairwise comparison






  