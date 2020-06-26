# Volcano plot to visualize GSEA results
# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html


# CREATED volcano_plot.R : 1/30/2020 by Cherry Gao
# MODIFIED for GSEA : 2/6/2020



#####  install  #####
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano')
#####################

# load package into R session
# library(EnhancedVolcano)
# library(stringr)
# library(dplyr)

##### START
rm(list = ls()) # clear environment
mother_dir = '/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/gsea/gsea_results/kegg_pathway_genesets_only'
setwd(mother_dir)
list_of_files = list.files() # list of folders in "gsea_results"
useful_cols = c(1,4:11) # manual; useful columns in GSEA result table 

# set thresholds for volcano plot
fdr_thresh = 0.25 # threshold on FDR q-value
nes_thresh = 1 # threshold on NES


for (z in 1:length(list_of_files)){
  # parse
  list_of_files[z] # print the folder to be processed
  setwd(paste0(mother_dir,'/',list_of_files[z]))
  neg = list.files(pattern = "^gsea_report_for_na_neg_(.*)xls$") # files that start with "gsea_report_for..." and ends with "xls"
  pos = list.files(pattern = "^gsea_report_for_na_pos_(.*)xls$") # files that start with "gsea_report_for..." and ends with "xls"
  pair_to_compare = unlist(strsplit(list_of_files[z],"\\."))[1] # split on '.' and extract the pair that is being compared
  
  # load GSEA result tables
  negative = read.table(file=neg, sep = "\t", header = TRUE) 
  positive = read.table(file=pos, sep = "\t", header = TRUE) 
  
  # get rid of useless columns
  negative = negative[,useful_cols]
  positive = positive[,useful_cols]
  
  # concatenate negative and positive
  data_concat = rbind(negative,positive)
  dim(positive)[1] + dim(negative)[1] == dim(data_concat)[1] # check number of rows
  
  # volcano plot
  pdf(file=paste(mother_dir,'/volcano_', pair_to_compare, '_fdrthresh_',fdr_thresh,'_nesthresh_',nes_thresh,'.pdf',sep="")) # to save
  
  print(EnhancedVolcano(data_concat,   # print needed to save pdf
                  lab = data_concat$NAME,       # dot labels
                  #                selectLab = c('647170520','647171017'),  # only label select points - text must also be contained in 'lab'
                  x = 'NES',          # data on x axis
                  y = 'FDR.q.val',                # data on y axis
                  xlim = c(-4, 4), 
                  ylim = c(-0.1, 5),
                  xlab = 'normalized enrichment score (NES)',
                  ylab = '-log10(FDR q-value)',
                  title = paste(pair_to_compare, '(FDR q-val=',fdr_thresh,'; NES=',nes_thresh,')'),
                  pCutoff = fdr_thresh,           # raw p-value cutoff threshold 
                  FCcutoff = nes_thresh,            # log2 fold change cutoff threshold
                  pointSize = 2,             # dot size
                  labSize = 1.5,              # text label size; 0 = no label
                  #                col=c('black','black','black','red'),       # in the order of NS -> log2FC -> p-value -> p-value + log2FC
                  legend=c('ns','NES','FDR q-val','FDR q-val + NES'),
                  legendVisible = FALSE,
                  colAlpha = 0.3,
                 drawConnectors = TRUE,   # connectors for labels
                 widthConnectors = 0.2,
                 colConnectors = 'grey30'))

    dev.off() # to save image
}





#################### OLD from volcano_plot.R

# extract the values of cutoff to inspect list of genes
# i_genes = which(abs(data_concat$foldchange) > fc_thresh & data_concat$padj < p_thresh) # desired index
# length(i_genes) # check length
# data_concat$geneProduct[i_genes] # check genes

# extract data and prep for GSEA
# select = list() # initialize
# select$KO = data_concat$KO[i_genes] # K numbers
# select$foldchange = data_concat$foldchange[i_genes] # fold change (for ranking)
# select = as.data.frame(select) # needed for 'mutate' function
#   select <- select %>% mutate(KO = str_extract(KO, "(?<=KO:)[0-9A-Za-z]*")) # extract K numbers only
#   select <- select %>% filter(!is.na(KO)) # get rid of NA (genes with no K assignment)
# head(select) # check
# length(select$KO) # check

# save_fn = paste0('gsea_input/test_rank.rnk')
# fn_save = paste0(cond1,'_vs_',cond2,'_select_pthresh_',p_thresh,'_fcthresh_',fc_thresh,'.rnk')
# write.table(as.data.frame(select), fn_save, append = FALSE, sep = "\t",quote=F, row.names = F, col.names = F)

