# Rank genes for GSEA (using UI)
# Create and save .rnk files for each pairwise comparison from DESeq2.


#----------SUMMARY----------
# Rank genes by log2FC - this was chosen because log2FC is biologically meaningful, and I have an alpha cutoff in DESeq so I'm kind of combining p-value and log2FC.
  ## Gene trimming steps: 
  # 1) Start with 5022 genes in DESeq output. 
  # 2) Then trim by padj cutoff (alpha = 0.05).
  # 3) Save rank table (column 1 = JGI gene numbers; column 2 = log2 fold change)
#---------------------------


#----------HISTORY----------
# UPDATED : 2/6/2020
#           added 10cont vs. 0cont pairwise comparison
# UPDATED : 3/17/2020 to v2
#           Instead of KO numbers, now using JGI gene ID numbers.
#           Automatically loop through each pairwise comparison.
#           No longer have the problem of multiple genes assigned.
#---------------------------

#-----INSTALL
library(dplyr)
library(stringr)



#-----SETUP
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load('pairwise_list.Rdata') # variable 'cond_pair_list'

# specify DESeq results folder to load
fdr_alpha = 0.01
folder = paste0('gsea/deseq_kegg_concat/')

# make a list of .Rdata files to load (DESeq pairwise comparisons), FULL directory
  # AND list of save names for the .rnk files (including full directory)
deseq_result_fn = matrix() # initialize
gsea_input_save_fn = matrix() # initialize
for (i in 1:length(cond_pair_list)) {
  deseq_result_fn[i] = paste0(folder, cond_pair_list[[i]][1],'_vs_',cond_pair_list[[i]][2],'_alpha_', fdr_alpha,'_KEGGconcat.Rdata')
  gsea_input_save_fn[i] = paste0('gsea/gsea_input/', cond_pair_list[[i]][1],'_vs_',cond_pair_list[[i]][2],'_alpha_', fdr_alpha,'.rnk')
}


# padj cutoff to apply to DESeq2 results before .rnk file creation
padj_cutoff = 0.01 

# create .rnk file for each pairwise comparison
for (z in 1:length(deseq_result_fn)) {
  
  # load data_concat - includes all genes including with high padj
  load(deseq_result_fn[z]) 
    # head(data_concat$foldchange)
    # head(data_concat$padj)
  
  # initialize
  rank = list()
  rank_trimmed = list()
  
  # extract variables
  rank = select(data_concat, seq_data_gene_id, foldchange, padj) # extract JGI gene ID#, log2FC, padj columns only
    # head(rank)
    # length(rank[,1])
  
  # change 'factor' to 'numeric'
    # foldchange and padj are a factor -- change to numeric
  rank$foldchange <- as.numeric(as.character(rank$foldchange)) # need to first convert to string, then convert foldchange to numeric 
    # print(as.numeric(paste(rank$foldchange))[1],digits=20) # at first, it looks like this factor -> numeric operation rounds the number but it actually doesn't (run this line to check)
  rank$padj <- as.numeric(as.character(rank$padj)) 
  
  # TRIM: apply padj cutoff & get rid of padj = NA (most likely small gene counts)
  nonsignificant_i = which(rank$padj>padj_cutoff | is.na(rank$padj)) # indexes of padj above cutoff
  rank_trimmed = rank[-c(nonsignificant_i),] # get rid of indexes with high pdaj
  
  # extract only JGI ID# and fold change as rank file, ready for GSEA
  rank_final = select(rank_trimmed, seq_data_gene_id, foldchange)
  
  # SAVE .rnk file for GSEA UI software
  write.table(as.data.frame(rank_final), gsea_input_save_fn[z], append = FALSE, sep = "\t", quote=F, row.names = F, col.names = F)
  
}

