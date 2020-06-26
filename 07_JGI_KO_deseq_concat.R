# Concatenate JGI gene ID, KO#, and DESeq2 results

# Feed in cleaned dataset generated from JGI_KO_dataset_clean.R.
# Create the variable, 'data_concat' which includes gene info from JGI database, concatenated with DESeq2 results (e.g. padj and fold change).
# Save 'data_concat' into folder '/gsea/deseq_kegg_concat'

#----------SUMMARY----------
# 1) load DeSeq expression results 
# 2) match JGI gene ID numbers and concatenate spreadsheets
# 3) save: resulting spreadsheet contains DeSeq results + Vcor genome data (incl. KO numbers) from JGI
#---------------------------


#----------HISTORY----------
# CREATED : 1/2020 by Cherry Gao
#           Note: the script works -- but not optimized (many counts are manually done)
# UPDATED : 3/16/2020
#           Instead of saving KO numbers, assign each JGI gene ID to a KEGG pathway to overcome multiple JGI gene IDs being assigned to a KO number.
#           Rename script from "deseq2_kegg_concat.R" to "JGI_KO_deseq_concat.R"
#---------------------------



# setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')

# load cleaned Vcor genome data downloaded from JGI (created from JGI_KO_dataset_clean.R)
load('dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # variable 'data_jgi'
load('pairwise_list.Rdata') # variable 'cond_pair_list'

# specify DESeq results folder to load
fdr_alpha = 0.01
deseq_result_folder = paste0('deseq/deseq2_results_alpha_',fdr_alpha,'_logfc_cutoff_none/')

# make a list of .Rdata files to load (DESeq pairwise comparisons)
deseq_result_fn = matrix() # initialize
for (i in 1:length(cond_pair_list)) {
  deseq_result_fn[i] = paste0(cond_pair_list[[i]][1],'_vs_',cond_pair_list[[i]][2],'_alpha_', fdr_alpha)
}


# extract first element of each list in data_jgi (i.e. only JGI gene ID numbers)
jgi_gene_list = unlist(lapply(data_jgi,'[[',1)); 

# for (z in 1:length(deseq_result_fn)) {
for (z in c(7:8)) {
    
  # read DESeq result file
  full_file_name_deseq = paste0(deseq_result_folder,deseq_result_fn[z],'.csv')
  deseq_result = readLines(full_file_name_deseq) # character object; refer to each line as line = raw_jgi[line_number]
  head(deseq_result) # look at table
  
  # deseq_result = read.csv(full_file_name_deseq)

  data_concat = matrix(nrow=length(deseq_result), ncol=34)  # list() # initialize
  
  
  for (j in 2:length(deseq_result)) {
    line = deseq_result[j] # just for test
    
    # handle non-header lines
    #  else {
    # line = unlist(strsplit(line, ',\t')) # read comma separated csv
    line = unlist(strsplit(gsub("\"","",line) , ",")) # get rid of extra string, then comma separate
    if (length(line) == 11) {
      line[10] =paste(line[10],line[11],sep=',') # paste back together
      line[11] = '' # make empty
    } # split too many times
    gene_oid = line[2] # JGI gene ID
    #   }
    
    # concatenate JGI genome data and expression data
    if (gene_oid %in% jgi_gene_list) {
      i = match(gene_oid,jgi_gene_list)
      if (gene_oid != data_jgi[[i]][1]){'expression data gene ID is not equal to JGI gene ID!'} # sanity check
      jgi_columns = data_jgi[[i]]  
      new_line = paste(c( line[1:2] , jgi_columns , line[3:10] ), sep = '\t') # concatenate headers
      # data_concat[[j]] = new_line # iteratively collect new lines
      # data_concat[3,] <- matrix(ncol=34, nrow=1)
      data_concat[j,] <- t(new_line)#rbind(data_concat,new_line) # iteratively collect new lines
    }
  } # end of cycling through line numbers
  
  # handle header -- first as the first row
  # if (j == 1) {
  # split the line into a list   #  line = line.rstrip().split('\t')
  line = deseq_result[1] 
  line = unlist(strsplit(line, ',')) # read comma separated csv
  line[1] = 'row_number' # rename; this is the row number in spreadsheet -- meaningless
  line[2] = 'seq_data_gene_id' # rename; originally 'rownames', change to 'seq_data_gene_id'
  line = unlist(strsplit(gsub("\"","",line) , ",")) # get rid of extra string, then comma separate
  
  header = paste(c( line[1:2] , columns , line[3:length(line)] ), sep = '\t') # concatenate headers
  data_concat[1,] = header 
  
  # header convert to names
  colnames(data_concat) = data_concat[1,] # make first row into colnames
  data_concat <- data_concat[-1,] # get rid of the first row (which are now headers)
  data_concat <- as.data.frame(data_concat) # convert to data frame -- now can refer to each variable with $
  
  # SAVE concatenated result to file
  save_fn = paste0('deseq/deseq_kegg_concat/',deseq_result_fn[z],'_KEGGconcat.txt')
  write.table(as.data.frame(data_concat), save_fn, append = FALSE, sep = "\t",quote=F, row.names = F, col.names = T)
  
  # SAVE to .RData (so that I don't have to read table)
  save_fn_r = paste0('deseq/deseq_kegg_concat/',deseq_result_fn[z],'_KEGGconcat.Rdata')
  save(data_concat,file = save_fn_r)
} # end of cycling through files
