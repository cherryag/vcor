# Concatenate t = 0 results

#---------- HISTORY ----------
# CREATED : 9/7/2020 by Cherry Gao
#-----------------------------

#----------SUMMARY----------
# 1) Load 10 muc vs. 0 muc and 60 muc vs. 0 muc results (BASELINE RESULTS)
# 2) Load 10 muc vs. 10 cont and 60 muc vs. 60 cont results (MAIN RESULTS)
# 3) Concatenate by Locus_tag (VIC numbers)
#---------------------------


library(reshape)
library(ggplot2)
library(scales)
library(readxl)


#---------- SET UP
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
  wrkdir = getwd()
data_dir = paste0(wrkdir,'/deseq/deseq_EEX_concat')
save_dir=paste0(wrkdir,'/deseq/deseq_0min_results_concat') # where to save new data sets
# dir.create(save_dir)

#---------- 1) Load BASELINE (vs. 0 min) results
# Load 10 muc vs. 0 muc and 60 muc vs. 0 muc results 

# 10 min base
load(file=paste0(data_dir,'/10muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))
base_data_concat_10 = data_concat
rm(data_concat)

# 60 min base
load(file=paste0(data_dir,'/60muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))
base_data_concat_60 = data_concat
rm(data_concat)

#---------- 2) Load MAIN (vs. time point matched) results
# Load 10 muc vs. 10 cont and 60 muc vs. 60 cont results 

# 10 min main
load(file=paste0(data_dir,'/10muc_vs_10cont_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))
main_data_concat_10 = data_concat
rm(data_concat)

# 60 min main
load(file=paste0(data_dir,'/60muc_vs_60cont_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))
main_data_concat_60 = data_concat
rm(data_concat)


#---------- 3) MATCH by EEX numbers
gene_idx_10 = match(main_data_concat_10$eex, base_data_concat_10$eex)
gene_idx_60 = match(main_data_concat_60$eex, base_data_concat_60$eex)


# check matching
i = 1000 # MANUAL
  base_data_concat_10$eex[gene_idx_10[i]] 
  main_data_concat_10$eex[i]

  base_data_concat_60$eex[gene_idx_60[i]] 
  main_data_concat_60$eex[i]

  
  
#---------- 4) Append EEX and save
  
# APPEND
main_data_concat_10$muc0_foldchange = base_data_concat_10$foldchange[gene_idx_10]
main_data_concat_10$muc0_padj = base_data_concat_10$padj[gene_idx_10]
main_data_concat_10$muc0_eex = base_data_concat_10$eex[gene_idx_10]

main_data_concat_60$muc0_foldchange = base_data_concat_60$foldchange[gene_idx_60]
main_data_concat_60$muc0_padj = base_data_concat_60$padj[gene_idx_60]
main_data_concat_60$muc0_eex = base_data_concat_60$eex[gene_idx_60]

# 
main_data_concat_10_appended = list()
main_data_concat_60_appended = list()
for (i in 1:2){
  
  # define data
  data_concat = list()
  if (i==1) {data_concat = main_data_concat_10
  } else if (i==2) {data_concat = main_data_concat_60}
    
    # loop through each gene
    for (j in 1:nrow(data_concat)) {
      
      if (data_concat$padj[j] == 'NA'){
        data_concat$main_base_results_match[j] = NA
        
      } else {
      
            # new label column: indicate significantly MAIN DE genes with "up" or "down"
            if (as.numeric(paste0(data_concat$padj[j])) < 0.01 && as.numeric(paste0(data_concat$foldchange[j])) > 0) {
                data_concat$significant_expression[j] = 'up'
            } else if (as.numeric(paste0(data_concat$padj[j])) < 0.01 && as.numeric(paste0(data_concat$foldchange[j])) < 0) {
                 data_concat$significant_expression[j] = 'down'
            } else {data_concat$significant_expression[j] = 'ns_main'}
            
            # new label column: indicate significantly BASELINE DE genes with "up" or "down"
            #if (as.numeric(paste0(data_concat$muc0_padj[j])) < 0.01 && as.numeric(paste0(data_concat$muc0_foldchange[j])) > 0) {
            if (as.numeric(paste0(data_concat$muc0_foldchange[j])) > 0) {
              data_concat$muc0_significant_expression[j] = 'up'
            #} else if (as.numeric(paste0(data_concat$muc0_padj[j])) < 0.01 && as.numeric(paste0(data_concat$muc0_foldchange[j])) < 0) {
            } else if (as.numeric(paste0(data_concat$muc0_foldchange[j])) < 0) {
              data_concat$muc0_significant_expression[j]= 'down'
            } else {data_concat$significant_expression[j] = 'ns_baseline'}
              
            
            # Do MAIN and BASELINE results match?
            if (data_concat$significant_expression[j] == data_concat$muc0_significant_expression[j]) {
              data_concat$main_base_results_match[j] = 'yes'
            } else if (data_concat$significant_expression[j] != data_concat$muc0_significant_expression[j]) {
              data_concat$main_base_results_match[j] = 'no'}
      }
  }

  # Redefine data
  if (i==1) {main_data_concat_10_appended = data_concat
  } else if (i==2) {main_data_concat_60_appended = data_concat}
}
      

# inspect genes that do not match MAIN and BASELINE results
sum(main_data_concat_10_appended$main_base_results_match == 'yes', na.rm=TRUE)
sum(main_data_concat_60_appended$main_base_results_match == 'yes', na.rm=TRUE)

# eliminated
num_elim_10 = sum(main_data_concat_10_appended$main_base_results_match == 'no', na.rm=TRUE)
num_elim_60 = sum(main_data_concat_60_appended$main_base_results_match == 'no', na.rm=TRUE)

# upregulated
num_up_10 = sum(main_data_concat_10_appended$significant_expression == 'up' & main_data_concat_10_appended$main_base_results_match == 'yes', na.rm=TRUE)
num_up_60 = sum(main_data_concat_60_appended$significant_expression == 'up' & main_data_concat_60_appended$main_base_results_match == 'yes', na.rm=TRUE)

# downregulated
num_down_10 = sum(main_data_concat_10_appended$significant_expression == 'down' & main_data_concat_10_appended$main_base_results_match == 'yes', na.rm=TRUE)
num_down_60 = sum(main_data_concat_60_appended$significant_expression == 'down' & main_data_concat_60_appended$main_base_results_match == 'yes', na.rm=TRUE)

# compute percentages  
num_up_10 / 5022 * 100
num_down_10 / 5022 * 100
num_up_60 / 5022 * 100
num_down_60 / 5022 * 100

# what I had before directionality match
sum(main_data_concat_10_appended$significant_expression == 'up')
sum(main_data_concat_10_appended$significant_expression == 'down')

sum(main_data_concat_60_appended$significant_expression == 'up')
sum(main_data_concat_60_appended$significant_expression == 'down')


#  SAVE as .Rdata and Excel (10 min DESeq2 results)
write.csv(main_data_concat_10_appended, file = paste0(save_dir,'/10muc_vs_10cont_alpha_0.01_KEGGconcat_EEXconcat_0minResConcat.csv'))
save(main_data_concat_10_appended, file = paste0(save_dir,'/10muc_vs_10cont_alpha_0.01_KEGGconcat_EEXconcat_0minResConcat.Rdata'))
  
# SAVE as .Rdata and Excel (60 min DEseq2 results)
write.csv(main_data_concat_60_appended, file = paste0(save_dir,'/60muc_vs_60cont_alpha_0.01_KEGGconcat_EEXconcat_0minResConcat.csv'))
save(main_data_concat_60_appended, file = paste0(save_dir,'/60muc_vs_60cont_alpha_0.01_KEGGconcat_EEXconcat_0minResConcat.Rdata'))