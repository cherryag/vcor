# Process manually curated virulence gene list and prepare for plotting

# This list is manually curated on an .xlsx through literature research.
# Should be run after each update on the Excel sheet, to convert to .Rdata file.

# CREATED : 3/23/2020 by Cherry Gao
# UPDATED : 4/25/2020
#           Process manually curated virulence genes after confirmation with NCBI accessions.

# library("readxl")

#---------- set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data



# convert data_jgi (maybe easier to work with)
data_jgi = do.call("rbind",data_jgi)
colnames(data_jgi) = unlist(columns) # append column names

# read in manually curated virulence genes
# v = read_excel('dataset/virulence_genes_curated/virulence_genes_vcor_trimmed.xlsx') # load dds
# v = read_excel('dataset/virulence_genes_curated/virulence_genes_vcor_trimmed_ncbi.xlsx') # load dds [ADDED 4/25/2020]
v = read_excel('dataset/virulence_genes_curated/vps_genes.xlsx') # load dds [ADDED 5/2/2020]
v = read_excel('dataset/virulence_genes_curated/secretion_system_genes.xlsx') # load dds [ADDED 5/2/2020]

jgi_id_and_name = cbind(v$`JGI #`,v$`gene name`)
# clean up list of JGI gene ID and gene names
jgi_id_and_name[,1]
index_to_eliminate = matrix() # initialize the list
  index_to_eliminate = which(is.na(jgi_id_and_name[,1]))
  index_to_eliminate = c(index_to_eliminate,  which(jgi_id_and_name[,1] == 'na'))
jgi_id_and_name = jgi_id_and_name[-index_to_eliminate,] # eliminated genes with no JGI gene ID assignment

# further cleanup list: separate 2 entries of JGI gene IDs into separate rows with repeated gene names
genes_for_timelapse_plot = list()
for (z in 1:length(jgi_id_and_name[,1])){
  
  ids = unlist(strsplit(jgi_id_and_name[z,1],", ")) # first column contains JGI gene ID (may be multiple)
  gene_name = jgi_id_and_name[z,2] # second column is gene name
  
  # if more than 1 JGI gene IDs, 
  if (length(ids)>1){ 
    # make a new row with gene name for each JGI gene ID
    for (i in 1:length(ids)) {     
      genes_for_timelapse_plot = rbind(genes_for_timelapse_plot, c(ids[i], gene_name))}
    # if only 1 JGI gene ID, just append the row as is
  } else if (length(ids) == 1) {genes_for_timelapse_plot = rbind(genes_for_timelapse_plot, c(ids,gene_name))}
}


# check result (manually curated virulence genes)
genes_for_timelapse_plot
jgi_id_and_name_manual = genes_for_timelapse_plot # rename for consistency with other sets of virulence genes

# SAVE
# save(jgi_id_and_name_manual, file='dataset/virulence_genes_curated/jgi_gene_id_virulence_manual_list_trimmed_ncbi.Rdata')
save(jgi_id_and_name_manual, file='dataset/virulence_genes_curated/jgi_gene_id_vps.Rdata')
