# Find JGI gene IDs of virulence genes by keywords

# in preparation for timelapse gene expression plots.

# CREATED : 3/23/2020 by Cherry Gao
# UPDATED : 5/13/2020 
#           Look for secretion system genes

### setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata')# load cleaned & concatenated JGI gene ID data

# convert data_jgi (maybe easier to work with)
data_jgi = do.call("rbind",data_jgi)
colnames(data_jgi) = unlist(columns) # append column names


# MANUALLY curated list (from literature search, not listed in virulence_genes.xlsx)
keywords_KO = c('K02453') # c('cholera toxin','vibriolysin') # keywords in KO only
keywords_product = c('zinc metalloprotease', 'MSHA', 'hemolysin') # keywords in product only
keywords_KO_COG_product = c('serine protease','NADH:ubiquinone') # keywords in KO, COG, and product 
keywords_KEGG_module = # keywords in KEGG module only [added 5/13/2020]


#----- keywords in KO only
jgi_id_and_name_KO = list()
for (i in 1:length(keywords_KO)){
  keyword = keywords_KO[i]
  col_index = which(columns == 'KO')

  # find
  index = c(grep(keyword, data_jgi[,col_index])) 
  index = unique(index)
  
  # store
  jgi_id_and_name_KO[[keyword]] = cbind(data_jgi[index,1], data_jgi[index,5]) # JGI gene IDs and gene product names
}

#----- keywords in product only
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

#-----keywords in KO, COG, and product 
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





#---- SAVE
save(jgi_id_and_name_KO,jgi_id_and_name_prod,jgi_id_and_name_KO_COG_prod,file='dataset/virulence_genes_curated/jgi_gene_id_virulence_keywords.Rdata')




# NEXT: 13_time_vs_gene_expression.R