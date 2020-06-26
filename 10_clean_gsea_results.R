# Clean up GSEA analysis results and save as .RData
# ONLY run once - to generate .RData

# CREATED : 3/18/2020 by Cherry Gao



############### SET UP ###############
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load('pairwise_list.Rdata')
deseq_alpha = 0.01
clean_results_savedir = paste0('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/gsea/gsea_results/fdr_',deseq_alpha)

# GSEA results folder
gsea_raw_results_dir = paste0('gsea/gsea_results/fdr_',deseq_alpha,'/raw_results')
setwd(gsea_raw_results_dir)
result_list = dir() # list of gsea results folders
result_list

# access each pairwise result folder, read table of results, and SAVE as .Rdata
for (z in 2:length(cond_pair_list)){ # skip 0m vs. 0c
  
  # initialize
  rm(neg_genesets_table,pos_genesets_table)
  
  # access pairwise comparison folder
  pairwise_compare = paste0(cond_pair_list[[z]][1],'_vs_',cond_pair_list[[z]][2])
  index_folder = pmatch(pairwise_compare, result_list)
  
  folder_name = result_list[index_folder]
  results_number = unlist(strsplit(folder_name,"\\."))[3] # long string of number, split with "."
  
  neg_genesets = paste0(folder_name,'/gsea_report_for_na_neg_',results_number,'.xls') # excel file
  pos_genesets = paste0(folder_name,'/gsea_report_for_na_pos_',results_number,'.xls') # excel file
  
  neg_genesets_table = read.table(file=neg_genesets, sep="\t", header=TRUE)
  neg_genesets_table = neg_genesets_table[,-c(2,3)] # get rid of unnecessary columns
  
  pos_genesets_table = read.table(file=pos_genesets, sep="\t", header=TRUE)
  pos_genesets_table = pos_genesets_table[,-c(2,3)] # get rid of unnecessary columns
  
  # save as Rdata
  save_name = paste0(clean_results_savedir,'/',pairwise_compare,'_gsea_results_clean.Rdata')
  save(deseq_alpha,pairwise_compare,results_number,neg_genesets_table,pos_genesets_table, file=save_name)
}
