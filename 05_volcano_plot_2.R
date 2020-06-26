# Volcano plot to visualize DESeq2 results (pairwise sample comparisons)
# If want to do subsequent GSEA on outliers on volcanos, need to go back to version 1 R script 
# where volcano plots were generated on KEGG-concatenated data. 
# (conclusion: GSEA on volcano outliers were not that interesting)
# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html


# CREATED : 1/30/2020 by Cherry Gao
# UPDATED : 3/14/2020 
#           Streamlined the code to take res_adj variable instead of post-KEGG assigned variable
#           For loop to iterate through pairwise comparisons


#####  install  #####
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano')
#####################

# load package into R session
library(EnhancedVolcano)
library(stringr)
library(dplyr)


##### START
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
alpha_change = 0.01 # FDR alpha used to generate DESeq2 results from pairwise comparisons
save_dir = paste0("deseq/deseq2_results_alpha_",alpha_change,"_logfc_cutoff_none/volcano_label")
# dir.create(save_dir)

# list of first condition (numerator; number on top)
cond_pair_list = list() # initialize
cond_pair_list[[1]] = c('0muc' , '0cont')
cond_pair_list[[2]] = c('10muc' , '10cont')
cond_pair_list[[3]] = c('60muc' , '60cont')
#cond_pair_list[[4]] = c('10cont' , '0cont')
#cond_pair_list[[5]] = c('60cont' , '0cont')
#cond_pair_list[[6]] = c('60cont' , '10cont')
#cond_pair_list[[7]] = c('10muc' , '0muc')
#cond_pair_list[[8]] = c('60muc' , '0muc')
#cond_pair_list[[9]] = c('60muc' , '10muc')

# set thresholds for volcano
p_thresh = 1e-15
fc_thresh = 0.5
w_fig = 10 # size of graphic region in inches (default value is 7)
h_fig = 10

# MANUAL list of genes to label
# gene_to_label = c('647172517','647171679','647173053','647173052','647169293','647169292') # Zn metalloprotease, VchA and B, ToxR and S
# gene_to_label = c('647172517') # Zn metalloprotease, VchA and B, ToxR and S
# gene_to_label = c('647172164','647172517','647173859') # metalloprotease genes in top 100 significantly DE genes at 60m vs. 60c
# gene_to_label = unlist(jgi_id_and_name_manual[,1]) # load in motility genes
gene_to_label = c('647169738', '647169740','647172193') # nqrs



for (j in 1:length(cond_pair_list)) { 
#  for (j in 1) { # to test
    
  # the 2 conditions to compare
  cond1 = cond_pair_list[[j]][1]
  cond2 = cond_pair_list[[j]][2]
  
  # load DESeq results
  fn = paste0('deseq/deseq2_results_alpha_',alpha_change,'_logfc_cutoff_none/',cond1,'_vs_',cond2,'_FDRalpha_',alpha_change,'.Rdata')
  load(fn) #load variable 'data_concat'
  
  # volcano plot
#  pdf(paste0(save_dir,'/volcano_', cond1, '_vs_', cond2, '_pthresh_',p_thresh,'_fcthresh_',fc_thresh,'.pdf'), w_fig,h_fig) # to save
  print(EnhancedVolcano(res_adj,   # print is needed for saving within for-loop
                  lab = rownames(res_adj), # res_adj$product, # paste(data_concat$Product_name,data_concat$jgi_gene_id, sep='_'),       # dot labels
  #                selectLab = c('647170520','647171017'),  # only label select points - text must also be contained in 'lab'
                  selectLab = gene_to_label,
                  x = 'log2FoldChange',          # data on x axis
                  y = 'padj',                # data on y axis
                  xlim = c(-10, 8), 
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  title = paste(cond1,'vs',cond2, '(p=',p_thresh,'; fc=',fc_thresh,')'),
                  pCutoff = p_thresh,           # raw p-value cutoff threshold 
                  FCcutoff = fc_thresh,            # log2 fold change cutoff threshold
                  pointSize = 1,             # dot size
                  labSize = 10,              # text label size; 0 = no label
                # col=c('black','black','black','red'),       # in the order of NS -> log2FC -> p-value -> p-value + log2FC
                  colAlpha = 0.3,
                  legend=c('ns','Log (base 2) FC','p-value','p-value + Log (base 2) FC'),
                  drawConnectors = TRUE, # connectors for labels -- not functioning as of 3/2020
                  widthConnectors = 0.5,
                  colConnectors = 'grey30', 
                  lengthConnectors = unit(0.01, 'npc') ) )       
#   dev.off() # to save image
  
  
  # extract list of genes that survived cutoffs
  i_genes = which(abs(res_adj$log2FoldChange) > fc_thresh & res_adj$padj < p_thresh) # desired index
  length(i_genes) # check length
  res_adj$product[i_genes] # check genes
  
  # save genes that survived cutoff (csv)
#  write.csv(as.data.frame(res_adj[i_genes,]),
#            file = paste0(save_dir,"/volcano_outliers_",cond1,"_vs_",cond2,"_padjthresh_",p_thresh,"_fcthresh_",fc_thresh,".csv"))

  
}