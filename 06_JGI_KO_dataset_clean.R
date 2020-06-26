# Assign JGI geneID to KEGG Orthology (KO) numbers

# Clean up the JGI database.
# Make the file, "jgi_gene_data_clean.Rdata" containing variables on interest (e.g. KO, gene length etc). 
# This script only has to be run once, since it's a database clean-up. 


# map JGI gene ID of Vcor to KO numbers (originally done by CS in python, contat.py)
# difference between 2016 and 2020: 
  # JGI database is slightly longer in 2020
# Downloaded Vcor (ATCC BAA-450) genes information from the JGI/IMG website (need to sign in) 
  # Click on "Export Gene Information" -> "Pages for download file 647000336.info.xls"
  # Manually change file extention to ".txt"
  # https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=GeneInfoPager&page=viewGeneInformation&taxon_oid=647000336


# 1) load Vcor genome data downloaded from JGI 
# 2) clean up data (reformat)
# 3) save workspace containing clean JGI data


# CREATED : 1/20/2020 by Cherry Gao (converted from CS's Python code (concat.py) to R)
# UPDATED : 3/16/2020 
#           just changed the script name fro "geneID2KO.R" to "JGI_KO_dataset_clean.R"
#           to do: also save as csv
#           to do: add header with each variable
# UPDATED : 3/23/2020
#           Fixed locus tag entry in jgi_data
#           Save jgi_data as csv file


# install.packages("tidyverse")



# set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/dataset')

# read JGI geneID table containing KO numbers (downloaded from internet)
file_name = 'jgi_downloaded/VcorBAA_genes_downloaded_from_JGI.txt'
raw_jgi = readLines(file_name) # character object; refer to each line as line = raw_jgi[line_number]
head(raw_jgi) # look at table
`%notin%` <- Negate(`%in%`) # create %notin% operator

# manually list column names to be extracted from Vcor genome data from JGI
columns <- list('jgi_gene_id', 'Locus_tag', 'Locus_type', 'Gene_symbol', 'Product_name', 
                'Scaffold', 'Coordinates', 'DNA_length', 'GC', 'COG_category', 
                'COG', 'pfam', 'TIGR', 'KEGG_module', 'KO', 
                'Metacyc', 'NCBI_accession', 'Protein_length', 'ITERM', 'EC', 
                'Transmembrane', 'Signal_peptide', 'Fused_gene', 'IMG')

# initialize 
data_jgi = list()
gene_id_list = c() # matrix with gene IDs only
line_count = 0 # line number, unique for each gene ID

for (j in 2:length(raw_jgi)) {    # start from row number 2 because of header
# for (j in 2:20) {    # for test
    line = raw_jgi[j] 
    
    if (line=="\t\t\t\t\t") {next} # skip empty line
  
    # split line into list (tab delimited)
    line = unlist(strsplit(line, '\t'))
      if (length(line) < 4) {next} # make sure the line has 4 entries (manually determined)
  
    # get information for this line
    gene_oid = line[1] # JGI gene ID is repeated over many lines
    locus_tag = line[2] # locus tag is repeated over many lines (The locus tag refers to the unique, systematic gene descriptor utilized by Genbank for genes from the V. coralliilyticus genome)
    ltype = line[3] # line type = type of information contained in this line
    other_cols = line[4:length(line)] # extract the other columns
    other_cols = other_cols[other_cols != ""] # discard empty
    linfo = paste(paste(other_cols,collapse=',')  ,sep=',') # join other column info with commas
    
    # add gene id and locus tag to our dictionary -- just once
    if (gene_oid %notin% gene_id_list) { 
          line_count = line_count + 1 # advance to new line number
          data_jgi[[line_count]] = matrix(nrow=1, ncol=length(columns)) # if container doesn't already have the gene ID, initialize with 'NA' for all columns
          data_jgi[[line_count]][1] = gene_oid # JGI gene ID is first column
          data_jgi[[line_count]][2] = locus_tag # locus tag is second clumn
          gene_id_list = c(gene_id_list,gene_oid) # iteratively add to gene ID list; separate vector just for the if statement 
    } 
  
    # for lines that perfectly match an entry in columns, become a column header
    if (ltype %in% columns) {
        i = match(ltype,columns) # index in columns that matches ltype
        data_jgi[[line_count]][i] = linfo
    # for lines that don't match column names
    } else {
        # handle KO number
        if (startsWith(ltype,'KO')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from 2 cols (in JGI data)
          ltype = 'KO'
        } else if (startsWith(ltype,'COG')) { # COG number
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'COG'
        } else if (startsWith(ltype,'pfam')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'pfam'
        } else if (startsWith(ltype,'TIGR')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'TIGR'
        } else if (startsWith(ltype,'ITERM')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'ITERM'
        } else if (startsWith(ltype,'EC')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'EC'
        } else if (startsWith(ltype,'IMG')) {
          linfo = paste(ltype,linfo, collapse=',', sep=',') # concatenate info from multiple cols (in JGI data)
          ltype = 'IMG'
        } else{'line entry does not match any of the column headers'} 
      
      # enter linfo at the right column (ltype)
      i = match(ltype,columns) 
      if (!is.na(data_jgi[[line_count]][i])) {   # if it's not NA, append linfo
         data_jgi[[line_count]][i] = paste(data_jgi[[line_count]][i],linfo,sep=',') # append linfo to existing entry
      } else {      # if it's NA, assign as linfo
        data_jgi[[line_count]][i] = linfo 
      }
    
    } # end of else
} # end of for every line loop

# convert to matrix for csv file
data_jgi_for_csv = do.call("rbind",data_jgi)
colnames(data_jgi_for_csv) = unlist(columns) # append column names

# SAVE as R file and .csv
save.image(file = "jgi_downloaded/jgi_gene_data_clean.Rdata")
write.csv(data_jgi_for_csv,file="jgi_downloaded/jgi_gene_data_clean.csv")


