# Process potential virulence genes from Kimes et al. 

# Extract locus tags of each gene set, and convert to JGI gene ID. 
# Source: Supplementary Table 1.  Spectral counting for all potential virulence factor proteins detected. (Kimes et al. 2012)


# CREATED : 3/23/2020 by Cherry Gao
# UPDATED : 3/25/2020 
#           Works with further edited Kimes Excel sheet (i.e. further break down QS genes)
# UPDATED : 4/18/2020
#           Include QS gene names (saved as _vs.Rdata file).
#           Extract proteases and save as .Rdata.

#---------- libraries
# install.packages("readxl")
library(readxl)
library(ggplot2)


#---------- Setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load data_jgi, cleaned & concatenated JGI gene ID data

# convert data_jgi (maybe easier to work with)
data_jgi = do.call("rbind",data_jgi)
colnames(data_jgi) = unlist(columns) # append column names

# read Kimes list of genes (cleaned .xlsx)
kimes = read_excel('dataset/virulence_genes_curated/Kimes_virulence_genes_cleaned.xlsx') # load dds

# list of gene sets in Kimes list [MANUAL] - must match exactly in the spreadsheet
kimes_geneset_list = c('Chemotaxis proteins',
                       'Methyl-accepting chemotaxis proteins (MCPs)',
                       'Flagellar proteins',
                       'Antibiotic resistance proteins',
                       'Hemolysins',
                       'Toxins',
                       'Proteases',
                       'General stress response proteins',
                       'Ribosomal proteins',
                       'Type 1 secretion proteins',
                       'Type 2 secretion proteins',
                       'Type four pilus proteins',
                       'Type III secretion proteins',
                       'Type VI secretion proteins',
                       'Quorum sensing (added)', # added 3/25/2020, formally 'Potential regulators'
                       'H-NS (added)', # added 3/25/2020, formally 'Potential regulators'
                       'Sigma factors (added)') # added 3/25/2020, formally 'Potential regulators'
                      #  'Potential regulators')

# loop through each gene set in Kimes list and extract locus tags, then convert to JGI gene list
jgi_gene_id_set_kimes = list() # initialize
# ALTERNATIVELY, don't loop through z and set z manually if specific gene group is desired. 
  z = 7 # proteases
for (z in 1:length(kimes_geneset_list)) {
  
  # determine start and end indexes of the gene set in Kimes list
  start = which(kimes$Category == kimes_geneset_list[z])
  end = which(kimes$Category == kimes_geneset_list[z+1]) # start of the NEXT category
  if (z==length(kimes_geneset_list)){end =length(kimes$Category)} # if end of Kimes list
  index_list = c(start:end) # list of indexes
  
  # list of locus tags + gene name list [added 4/18/2020]
  locus_tag_list = kimes[index_list,2] 
  gene_name_list = kimes[index_list,4]
  
  # cleanup list of locus tags
  index_to_eliminate = matrix() # iniitalize
  index_to_eliminate = which(locus_tag_list == 'N/A*') 
  index_to_eliminate = append(index_to_eliminate, which(is.na(locus_tag_list)) )
  
  # if list of index to eliminate is NOT empty, eliminate
  if (!length(index_to_eliminate)==0){ 
    locus_tag_list = locus_tag_list[-index_to_eliminate,]
    locus_tag_list = unlist(locus_tag_list) 
    # locus_tag_list # view the list
    
    gene_name_list = gene_name_list[-index_to_eliminate,] # [added 4/18/2020]
    gene_name_list = unlist(gene_name_list)
  }
  
  # match locus tag with JGI gene ID and save in permanent container
  jgi_index = match(locus_tag_list, data_jgi[,2])
  jgi_gene_id_set_kimes[[z]] = list() # initialize
  # jgi_gene_id_set_kimes[[z]] = as.numeric(data_jgi[jgi_index,1]) # save in permanent container
  jgi_gene_id_set_kimes[[z]] = cbind(data_jgi[jgi_index,1], gene_name_list) # [added 4/18/2020]
    rownames(jgi_gene_id_set_kimes[[z]]) <- c()
    colnames(jgi_gene_id_set_kimes[[z]]) <- c()
  
}

# assign gene set names ("z" only if not loop through)
temp[[kimes_geneset_list[[z]]]] = jgi_gene_id_set_kimes[[z]] # temp = dummy variable
jgi_gene_id_set_kimes <- temp # assign dummy variable to final data container  
# names(jgi_gene_id_set_kimes) = kimes_geneset_list # if loop through

# SAVE Kimes gene set JGI gene IDs
# save(jgi_gene_id_set_kimes,kimes_geneset_list,file='dataset/virulence_genes_curated/jgi_gene_id_sets_kimes_qs_v2.RData')
save(jgi_gene_id_set_kimes,kimes_geneset_list,file='dataset/virulence_genes_curated/jgi_gene_id_sets_kimes_protease.RData')


#----- NEXT : 12_bar_plot_gene_sets_by_jgi_id.R or 14_time_vs_log2fc.R

