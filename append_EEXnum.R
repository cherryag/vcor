# Append EEX protein numbers to data sets

#---------- HISTORY ----------
# CREATED : 7/16/2020 by Cherry Gao
#-----------------------------

#----------SUMMARY----------
# 1) Create EEX lookup table (saved as EEX_VIC_lookup.Rdata)
# 2) APPEND EEX to JGI gene data (whole genome)
# 3) APPEND EEX to DESeq2 results (for publishing)
#---------------------------


library(reshape)
library(ggplot2)
library(scales)
library(readxl)


#---------- SET UP
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
  wrkdir = getwd()
save_dir=paste0(wrkdir,'/dataset/EEX_appended') # where to save new data sets
# dir.create(save_dir)

#---------- 1) Create EEX lookup table (saved as EEX_VIC_lookup.Rdata)
# load dataset downloaded from NCBI on 5/14/2020, contains VIC numbers and EEX numbers for each gene
eex_dataset_dir = 'dataset/GenBank Vcor genome/ncbi-genomes-2020-05-14/GCA_000176135.1_ASM17613v1_feature_table.xlsx'
v = read_excel(eex_dataset_dir) 

# find indexes that contain EEX numbers (not NaN)
id_data = which(!is.na(v$product_accession)) 
  # check
  v$product_accession[id_data] # should contain all EEX
  length(id_data) # should be 5022 (# Vcor genes)
  
# construct look-up table for EEX numbers using VIC
eex_lookup = list() # initialize
eex_lookup$eex = v$product_accession[id_data] 
eex_lookup$vic = v$locus_tag[id_data]
eex_lookup$length_aa = v$product_length[id_data]
eex_lookup$name = v$name[id_data] 

# SAVE lookup table as .Rdata
# save(eex_lookup, file = paste0(save_dir,'/EEX_VIC_lookup.Rdata'))



#---------- 2) APPEND EEX to JGI gene data (whole genome)

# LOAD 
rm(v) # clean variable
load('dataset/jgi_downloaded/jgi_gene_data_clean.Rdata')

jgi_dataset_dir = 'dataset/jgi_downloaded/jgi_gene_data_clean.csv'
v = read.csv(jgi_dataset_dir,header=TRUE) 

# MATCH by VIC number (gene_idx = index in eex_lookup)
gene_idx = match(v$Locus_tag, eex_lookup$vic)

# check matching
i = 5000 # MANUAL
eex_lookup$vic[gene_idx[i]] 
  v$Locus_tag[i]

# APPEND
v$eex = eex_lookup$eex[gene_idx]

# CLEAN
v$X <- NULL # get rid of the extra column
  
# SAVE as .Rdata and Excel
write.csv(v, file = paste0(save_dir,'/jgi_gene_data_clean_eex_appended.csv'))



#---------- 3) APPEND EEX to DESeq2 results (for publishing)
# LOAD EEX lookup table
load(paste0(save_dir,'/EEX_VIC_lookup.Rdata'))

# load DESeq2 results, 10 min
#load(file='deseq/deseq_kegg_concat/10muc_vs_10cont_alpha_0.01_KEGGconcat.Rdata')
load(file='deseq/deseq_kegg_concat/10muc_vs_0muc_alpha_0.01_KEGGconcat.Rdata')

data_concat_10 = data_concat
rm(data_concat)

# load DESeq2 results, 60 min
#load(file='deseq/deseq_kegg_concat/60muc_vs_60cont_alpha_0.01_KEGGconcat.Rdata')
load(file='deseq/deseq_kegg_concat/60muc_vs_0muc_alpha_0.01_KEGGconcat.Rdata')

data_concat_60 = data_concat
rm(data_concat)

# MATCH
gene_idx_10 = match(data_concat_10$Locus_tag, eex_lookup$vic)
gene_idx_60 = match(data_concat_60$Locus_tag, eex_lookup$vic)

#gene_idx_60 = match(as.numeric(unlist(gene_list[,1])), data_concat_60$seq_data_gene_id)

# check matching
i = 5000 # MANUAL
  eex_lookup$vic[gene_idx_10[i]] 
  data_concat_10$Locus_tag[i]
  
  eex_lookup$vic[gene_idx_60[i]] 
  data_concat_60$Locus_tag[i]

# APPEND
data_concat_10$eex = eex_lookup$eex[gene_idx_10]
data_concat_60$eex = eex_lookup$eex[gene_idx_60]

#  SAVE as .Rdata and Excel (10 min DESeq2 results)
data_concat = data_concat_10 # convert back
write.csv(data_concat, file = paste0(save_dir,'/10muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.csv'))
save(data_concat, file = paste0(save_dir,'/10muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))

# SAVE as .Rdata and Excel (60 min DEseq2 results)
rm(data_concat)
data_concat = data_concat_60 # convert back
write.csv(data_concat, file = paste0(save_dir,'/60muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.csv'))
save(data_concat, file = paste0(save_dir,'/60muc_vs_0muc_alpha_0.01_KEGGconcat_EEXconcat.Rdata'))

