# Create top 50 genes table - DESeq2 results



# CREATED : 3/18/2020 by Cherry Gao
#           Coded up 2 methods of sorting: 
#             Method 1: Sort all subsampled (100) genes according to padj, then sort by positive or negative log2FC (1 table)
#             Method 2: Sort by padj, take the top genes (100), then sort by log2FC. 
#                       Create 2 tables: 40 genes that are positive and 40 genes that are negative.
#                       Decided against method 2 because it seems biased to take the subset of subset. 




#------- setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load('pairwise_list.Rdata')
deseq_alpha = 0.01
subset_gene_length = 100 # how many genes to sort by padj first
final_gene_length = 40 # how many genes to keep in each table
col_delete =   c(1,2,4,5,8,9,10,11,13,14,15:26,27,29,30,31,33,34) # index of columns to delete in the trimmed version of the table
save_fn = 'deseq/top_de_gene_table/' # directory to save tables

# full directory of 'data_concat' files to load
dir_n = 'deseq/deseq_kegg_concat/'
full_fn = matrix() # initialize
pairwise_vs = matrix() # initialize
for (z in 1:length(cond_pair_list)){
  rdata_fn = paste0(cond_pair_list[[z]][1],'_vs_', cond_pair_list[[z]][2],'_alpha_',deseq_alpha,'_KEGGconcat.Rdata')
  pairwise_vs[z] = paste0(cond_pair_list[[z]][1],'_vs_', cond_pair_list[[z]][2]) # for file name only
  full_fn[z] = paste0(dir_n,rdata_fn)
}


#------ METHOD 1: sort all subsampled (100) genes according to padj, then sort by positive or negative log2FC (1 table)
# iterate through each pairwise comparison and save top differentially expressed genes as tables
for (z in 1:length(cond_pair_list)){

  data_concat_trim = list() # initialize
  load(file=full_fn[z]) # load data_concat
  
  # data_concat housekeeping
  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
  data_concat$padj = as.numeric(paste(data_concat$padj))
  
  # data_concat_ordered_padj = list() # initialize
    # head(data_concat,10)
    # head(data_concat[order(data_concat$padj),],10)
  
  # ORDER 1: by padj
  data_concat_ordered_padj = data_concat[order(data_concat$padj),] # order by padj
  data_concat_subset = data_concat_ordered_padj[c(1:subset_gene_length),] # take the top subset of genes
    # head(data_concat_subset)
  
  
  # ORDER 2: order the subset of genes by log2 fold change
  data_concat_subset = data_concat_subset[order(-data_concat_subset$foldchange),] # order by fold change (positive to negative)
    # head(data_concat_subset)
  
  # TRIM: reorder columns (like final version of table)
  data_concat_trim$jgi_gene_id = data_concat_subset$jgi_gene_id
  data_concat_trim$jgi_gene_id = data_concat_subset$jgi_gene_id
  # delete redundant columns
  data_concat_trim = data_concat_subset[-col_delete]
  
  # SAVE tables
  write.table(data_concat_subset, file=paste0(save_fn,pairwise_vs[z],'_top_',subset_gene_length,'_genes_padj.csv'), sep=',', col.names=T, row.names = F)
  write.table(data_concat_trim, file=paste0(save_fn,pairwise_vs[z],'_top_',subset_gene_length,'_genes_padj_trim.csv'), sep=',', col.names=T, row.names = F)
}








#------ METHOD 2: Sort by padj, take the top genes (100), then sort by log2FC. Create 2 tables: 40 genes that are positive and 40 genes that are negative.
  # decided against this method because it seems biased to take the subset of subset. 

# iterate through each pairwise comparison and save top differentially expressed genes as tables
for (z in 1:legnth(cond_pair_list)){

  neg_fc_table = list() # initialize
  pos_fc_table = list() # initialize
  load(file=full_fn[z]) # load data_concat

  # data_concat housekeeping
  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
  data_concat$padj = as.numeric(paste(data_concat$padj))

  # data_concat_ordered_padj = list() # initialize
    # head(data_concat,10)
    # head(data_concat[order(data_concat$padj),],10)

  # ORDER 1: by padj
  data_concat_ordered_padj = data_concat[order(data_concat$padj),] # order by padj
  data_concat_subset = data_concat_ordered_padj[c(1:subset_gene_length),] # take the top subset of genes
    # head(data_concat_subset)


  # ORDER 2: by log2 fold change
  data_concat_subset = data_concat_subset[order(data_concat_subset$foldchange),] # order by fold change
    # head(data_concat_subset)


  # NEGATIVE DE genes table
  if (sum(!data_concat_subset$foldchange[1:final_gene_length] < 0) == 0){ # if none of the genes are NOT negative
      neg_fc_table = data_concat_subset[c(1:final_gene_length),]
  } else {print(paste(pairwise_vs[z],': need to increase the number of genes sampled - not all positive'))}

  # POSITIVE DE genes table (re-order)
  data_concat_subset = data_concat_subset[order(-data_concat_subset$foldchange),] # descending order, positive fold changes
    # head(data_concat_subset) 
  
  if (sum(!data_concat_subset$foldchange[1:final_gene_length] > 0) == 0){ # if none of the genes are NOT positive
      pos_fc_table = data_concat_subset[c(1:final_gene_length),]
  } else {print(paste(pairwise_vs[z],': need to increase the number of genes sampled - not all negative'))}

  # SAVE tables
  write.table(pos_fc_table, file=paste0(pairwise_vs[z],'_top_positive.csv'), sep=',', col.names=T, row.names = F)
  write.table(neg_fc_table, file=paste0(pairwise_vs[z],'_top_negative.csv'), sep=',', col.names=T, row.names = F)
}