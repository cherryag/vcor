# Validation ratio for validation of genes that DE in mucus
# To address Janelle's comments

#---------- HISTORY ----------
# CREATED : 7/17/2020 by Cherry Gao
#-----------------------------

#----------SUMMARY----------
# 1) load DESeq2 results for appropriate pairwise time comparisons, for mucus and control
# 2) calculate validation ratio (mucus / control)
#---------------------------

library(ggplot2)

#---------- SET UP
rm(list = ls()) # clear environment
wrkdir = '/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal'
setwd(wrkdir)
save_dir=paste0(wrkdir,'/deseq/time_de') # where to save timelapse plots
# dir.create(save_dir)


#---------- LOAD a pair of DESeq2 results for ratio calculation
tf = '10'
ti = '0'

load(paste0('deseq/deseq2_results_alpha_0.01_logfc_cutoff_none/',tf,'cont_vs_',ti,'cont_FDRalpha_0.01.Rdata'))
degenelist_c = degenelist
rm(degenelist)

load(paste0('deseq/deseq2_results_alpha_0.01_logfc_cutoff_none/',tf,'muc_vs_',ti,'muc_FDRalpha_0.01.Rdata'))
degenelist_m = degenelist
rm(degenelist)


# match genes
gene_idx_m = match(degenelist_c$rownames, degenelist_m$rownames)

# check match
i = 1000
degenelist_c$rownames[i] 
degenelist_m$rownames[gene_idx_m[i]]

# get id-matched log2 fold change of genes
m_l2 = as.numeric(paste(degenelist_m$foldchange[gene_idx_m]))
c_l2 = as.numeric(paste(degenelist_c$foldchange))

# convert log2 fold changes to fold changes
m = 2^m_l2
c = 2^c_l2

# calculate ratio
ratio <- m/c
#ratio <- na.omit(ratio)
#ratio <- ratio[which(is.finite(ratio))]
ratio_log10 = log10(ratio) # log10 transform, to make the histogram work (otherwise too long tail)

# histogram of ratios
hist(ratio_log10, breaks=30,
     main=paste0('mucus/control fold-change ratios \n T',tf,' vs T',ti, ' DESeq2 within mucus OR ctrl'),
     xlab='log10(muc / ctrl fold change)',
     ylab='counts')

# look at the list of genes
v = list()
v$jgi = degenelist_c$rownames
v$prod = degenelist_c$geneProduct
v$mcratio = ratio

write.csv(v,file=paste0(save_dir,'/mc_fc_ratio_',tf,'_vs_',ti,'.csv'))
