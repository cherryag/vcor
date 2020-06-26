# Determine number of differentially expressed genes that overlap between 10 and 60 min
# Used to generate Supplementary Table 3.

# CREATED : 8/7/2020 by Cherry Gao


#---------- set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/deseq/deseq2_results_alpha_0.01_logfc_cutoff_none')

# COPIED from 04_bar_plot_DE_genecount_alpha_0.01.R
up0 = 0 # 0M vs 0C
down0 = 0
up10 =  1379 # 10M vs 10C
down10 = 1326
up60 =  1159 # 60M vs 60C
down60 = 1076


#------------------------- load 10 min vs. 10 min results
load(file='10muc_vs_10cont_FDRalpha_0.01.Rdata') # DESeq results
degenelist_10 = degenelist
rm(degenelist)
  # all significantly DE genes
  id10 = which(as.numeric(paste0(degenelist_10$padj))<0.01)
  jgi10 = as.numeric(paste0(degenelist_10$rownames[id10])) # index, then JGI gene ID of padj > 0.01
  length(jgi10) == (up10+down10)
  # signifcantly upregulated genes
  id10_up = which(as.numeric(paste0(degenelist_10$padj))<0.01 & as.numeric(paste0(degenelist_10$foldchange))>0)
  jgi10_up = as.numeric(paste0(degenelist_10$rownames[id10_up])) # index, then JGI gene ID of padj > 0.01 & positive log2 fold change
  length(jgi10_up) == up10
  # significantly downregulated genes
  id10_down = which(as.numeric(paste0(degenelist_10$padj))<0.01 & as.numeric(paste0(degenelist_10$foldchange))<0)
  jgi10_down = as.numeric(paste0(degenelist_10$rownames[id10_down])) # index, then JGI gene ID of padj > 0.01 & negative log2 fold change
  length(jgi10_down) == down10

#-------------------------- load 60 min vs. 60 min results
load(file='60muc_vs_60cont_FDRalpha_0.01.Rdata') # DESeq results
degenelist_60 = degenelist
rm(degenelist)
  # all significantly DE genes
  id60 = which(as.numeric(paste0(degenelist_60$padj))<0.01)
  jgi60 = as.numeric(paste0(degenelist_60$rownames[id60])) # index, then JGI gene ID of padj > 0.01
  length(jgi60) == (up60+down60)
  # signifcantly upregulated genes
  id60_up = which(as.numeric(paste0(degenelist_60$padj))<0.01 & as.numeric(paste0(degenelist_60$foldchange))>0)
  jgi60_up = as.numeric(paste0(degenelist_60$rownames[id60_up])) # index, then JGI gene ID of padj > 0.01 & positive log2 fold change
  length(jgi60_up) == up60
  # significantly downregulated genes
  id60_down = which(as.numeric(paste0(degenelist_60$padj))<0.01 & as.numeric(paste0(degenelist_60$foldchange))<0)
  jgi60_down = as.numeric(paste0(degenelist_60$rownames[id60_down])) # index, then JGI gene ID of padj > 0.01 & negative log2 fold change
  length(jgi60_down) == down60


# Q1: How many genes are UPREGULATED at 10 min and UPREGULATED at 60 min? 
up_10_up_60 = length(intersect(jgi10_up,jgi60_up))
up_10_up_60/5020 * 100
  
# Q2: How many genes are DOWNREGULATED at 10 min and DOWNREGULATED at 60 min? 
down_10_down_60 = length(intersect(jgi10_down,jgi60_down))
down_10_down_60/5020 * 100

# Q3: How many genes are UPREGULATED at 10 min and DOWNREGULATED at 60 min? 
up_10_down_60 = length(intersect(jgi10_up,jgi60_down))
up_10_down_60/5020 * 100

# Q4: How many genes are DOWNREGULATED at 10 min and UPREGULATED at 60 min? 
down_10_up_60 = length(intersect(jgi10_down,jgi60_up))
down_10_up_60/5020 * 100

# [optional check] genes that switch between significant up and down between 10 and 60 min
genes_that_switch = setdiff(intersect(jgi10,jgi60),  c(intersect(jgi10_down,jgi60_down),intersect(jgi10_up,jgi60_up)))

# Q5: How many DE genes are shared between 10 and 60 min? 
all_i = up_10_up_60 + down_10_down_60 + up_10_down_60 + down_10_up_60
all_i == length(intersect(jgi10,jgi60))
all_i/5020 * 100








