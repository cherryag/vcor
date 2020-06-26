# Find gene information by keywords
# Modified from virulence_genes_keywords.R

# in preparation for timelapse gene expression plots.

# CREATED : 3/23/2020 by Cherry Gao (virulence_genes_keywords.R)
# UPDATED : 5/13/2020 to extract_genes_keywords.R
#           Look for secretion system genes
# DUPLICATED in 14_time_vs_log2fc_03.R for direct Excel to timelapse plot generation [5/28/2020].

# libraries
library("readxl")

### setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata')# load cleaned & concatenated JGI gene ID data

# convert data_jgi (maybe easier to work with)
data_jgi = do.call("rbind",data_jgi)
colnames(data_jgi) = unlist(columns) # append column names

# read in Excel sheet containing keywords
gene_list_of_interest = 'flagella' # MANUAL
# v = read_excel('dataset/virulence_genes_curated/secretion_system_genes.xlsx') # load dds [ADDED 5/2/2020]
v = read_excel(paste0('dataset/virulence_genes_curated/gene_list_',gene_list_of_interest,'.xlsx')) # load dds [ADDED 5/2/2020]

#-------------- (1) KEYWORDS = KO number or words -------------------------
# extract keywords
keywords = v$`KO number`

# clean up list of keywords
index_to_eliminate = matrix() # initialize the list
index_to_eliminate = which(is.na(keywords))
index_to_eliminate = c(index_to_eliminate,  which(keywords == '-'))
if(length(index_to_eliminate)>0){ keywords = keywords[-index_to_eliminate] } # eliminated 
keywords = unique(keywords) # no repeat

# what gene info to extract?
columns # print all variables
extract_cols = c(1:length(columns))

#-------- MANUALLY curated list (from literature search, not listed in virulence_genes.xlsx)
# keywords_KO = c('cholera toxin','vibriolysin') # keywords in KO only
# keywords_product = c('zinc metalloprotease', 'MSHA', 'hemolysin') # keywords in product only
# keywords_KO_COG_product = c('serine protease','NADH:ubiquinone') # keywords in KO, COG, and product 
# keywords_KEGG_module = # keywords in KEGG module only [added 5/13/2020]
#----------------------------------

#----- keywords in KO only
gene_info_extract = list()
for (i in 1:length(keywords)){
  key = keywords[i]
  col_index = which(columns == 'KO')
  
  # clean up key (if necessary)
  key = substring(key,1,6) # first digits = K numbers

  # find
  index = c(grep(key, data_jgi[,col_index])) 
  index = unique(index)
  
  # store
  #  gene_info_extract[[key]] = cbind(data_jgi[index,1], data_jgi[index,5]) # extract subset of gene info
  #gene_info_extract[[key]] = as.data.frame(data_jgi[index,extract_cols])
#   gene_info_extract_2[i,] = as.list(data_jgi[index,extract_cols])
  gene_info_extract = rbind(gene_info_extract, data_jgi[index,extract_cols])
}


#---- SAVE
write.table(gene_info_extract, file = paste0('gene_list_',gene_list_of_interest,'_extract.txt'), sep = "\t",row.names = FALSE)

#-------------- END of (1) KEYWORDS = KO number or words -------------------------



#----- keywords in product only [not updated]
jgi_id_and_name_prod = list()
for (i in 1:length(keywords_product)){
  keyword = keywords_product[i]
  col_index = which(columns == 'Product_name')
  
  # find
  index = c(grep(keyword, data_jgi[,col_index])) 
  index = unique(index)
  
  # store
  jgi_id_and_name_prod[[keyword]] = cbind(data_jgi[index,1], data_jgi[index,5]) # JGI gene IDs and gene product names
}

#-----keywords in KO, COG, and product [not updated]
jgi_id_and_name_KO_COG_prod = list()
for (i in 1:length(keywords_KO_COG_product)){
  
  keyword = keywords_KO_COG_product[i]
  
  p_index = which(columns == 'Product_name')
  ko_index = which(columns == 'KO')
  cog_index = which(columns == 'COG')
  
  # find
  index = c(grep(keyword, data_jgi[,p_index]), 
            grep(keyword, data_jgi[,ko_index]),
            grep(keyword, data_jgi[,cog_index]))
  
  # store
  jgi_id_and_name_KO_COG_prod[[keyword]] = cbind(data_jgi[index,1], data_jgi[index,5]) # JGI gene IDs and gene product names
  jgi_id_and_name_KO_COG_prod[[keyword]] = unique(jgi_id_and_name_KO_COG_prod[[keyword]])
}







