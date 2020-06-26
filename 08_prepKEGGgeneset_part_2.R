# Create geneset database defined by KEGG PATHWAY 
# PART 2 ONLY (clean up KEGG-downloaded database)

# Categorize genes (= KEGG Orthology (KO) numbers) for Vcor OCN014.
# Many of the gene counts for METHODS section is in this script.
# Save .gmt file for GSEA.


#--------------------SOURCE OF KEGG DATASET--------------------
# Downloaded as 'htext' (contains BRITE and PATHWAY) at https://www.genome.jp/kegg-bin/get_htext?vct00001
# Manually changed '.keg' file extension to '.txt' to enable loading into R. 
# NOTE: json file (and manually changing to .txt) gives KEGG PATHWAY but different format than 2016 file.
# NOTE: slight differences between 2016 and 2020 downloaded KEGG files
# 2016 downloaded KEGG file: 'VcorOCN014_KEGGpathways_original.unix.txt' also contained a mixture of BRITE and PATHWAY geneset categories.
# HOWEVER, in 2016, all BRITE categories were empty (didn't contain any genes, i.e. KO numbers) probably because BRITE was a relatively new database back then.
# Between 2016 and 2020, same KO numbers and associated gene names.
# Parsed geneset files: 2016: "PATH" only (because "BR" categories were empty); whereas 2020: mixture of "PATH" and BR" 
# Wikipedia: "Another database that supplements KEGG PATHWAY is the KEGG BRITE database. 
#             It is an ontology database containing hierarchical classifications of various entities including genes, 
#             proteins, organisms, diseases, drugs, and chemical compounds. While KEGG PATHWAY is limited to molecular 
#             interactions and reactions of these entities, KEGG BRITE incorporates many different types of relationships including:
#             Genes and Proteins; Compounds and Reactions; Drugs; Diseases; Organisms and Cells"
# KEGG BRITE database contains KEGG Orthology (KO), the reaction classification system for biochemical reactions, and other classifications for compounds and drugs.
# JV59_##### in downloaded KEGG file refer to genes in the Vcor genome.
#---------------------------------------------------------------

#----------SUMMARY----------
# PART 1: 
# 1) load KEGG file containing KO numbers (downloaded from KEGG website)
# 2) parse KEGG file: PATHWAY category numbers, description, and KO numbers
# 3) delete empty PATHWAY category numbers
# 4) save to file

# PART 2: 
# 1) replace KO# in the dataset created in Part 1 with JGI gene ID numbers
# 2) save to file (.gmt)
#---------------------------

#----------HISTORY----------
# CREATED v1 : 1/20/2020 by Cherry Gao
#              DOWNLOADED KEGG databases : 1/20/2020
# UPDATED to v2 : 2/6/2020  get rid of KEGG BRITE, and only use KEGG PATHWAY
#                 KEGG BRITE is too general -- even losely connected genes are grouped together. 
#                 We want the highest quality pathway groupings, i.e. KEGG PATHWAY.
#                 Also got rid of repeated KO numbers -- multiple KO numbers (with different Vcor gene IDs in downloaded KEGG list) may be assigned to a pathway.
# UPDATED v2 : 3/17/2020 
#              Created PART 2, where KO# are replaced by JGI gene ID # for subsequent GSEA. 
#              PART 2 can be stand-alone -- i.e. don't need to run Part 1 beforehand. 
#---------------------------



#----------PART 2: assign JGI gene ID numbers to KEGG PATHWAYs [ADDED 3/17/2020]----------
# Goal: create a spreadsheet with KEGG PATHWAY followed by columns of JGI gene ID#s
# e.g.: KEGG PATHWAY; JGI ID#; JGI ID#; JGI ID#...
# save to .gmt, file format for GSEA.


# set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')

# load(paste0('gsea/deseq_kegg_concat/',deseq_result_fn[z],'_KEGGconcat.Rdata'))
load('dataset/kegg_downloaded/kegg_pathways.Rdata')
load('dataset/jgi_downloaded/jgi_gene_data_clean.Rdata')

# indexes of variables of interest in 'data_jgi' (according to 06_JGI_KO_dataset_clean.R)
jgi_num_index = 1
gene_name_index = 5 # 'product name'
ko_num_index = 15 # 'KO'  


# initialize data containers
kegg_list_jgi = data.frame()
kegg_list_jgi = kegg_list # at first, append JGI gene ID# to KO#
count_pathways_match = matrix() # count the number of pathways to which a gene is assigned
no_ko_gene = matrix(ncol=2) # count the number of genes without KO # assignment
colnames(no_ko_gene) = c('gene name','JGI ID')
no_pathway_match_gene = matrix(ncol=3) # genes with no match to any pathways
colnames(no_pathway_match_gene) = c('gene name','JGI ID', 'KO')


# for each gene, iterate through each KEGG PATHWAY and append to PATHWAY when there is a match between KO# and JGI#
for (z in 1:length(data_jgi)) { # for all genes
  
  # extract JGI number and KO number of the gene
  jgi_num = data_jgi[[z]][jgi_num_index] # data_concat$jgi_gene_id[z]
  k_num = data_jgi[[z]][ko_num_index] # data_concat$KO[z]
  gene_name = data_jgi[[z]][gene_name_index] # data_concat$Product_name[z]
  
  # extract the KO number only (operations do nothing for 'NA')
  a = unlist(strsplit( paste(k_num), ','))[1]
  k_num = gsub('KO:*', '', a)  
  
  # if no assignment of KO number to JGI ID exists, keep track, and move onto the next gene
  if (k_num == 'NA') {
    #     no_ko_gene = append(no_ko_gene, gene_name)
    no_ko_gene = rbind(no_ko_gene, c(gene_name, jgi_num))
    count_pathways_match[z] = NA # assign NA to count_pathways_match[z]  (to not miss the last genes)
    next} 
  
  # initialize count
  count_pathways_match_temp = matrix()
  
  # iterate through each KEGG PATHWAY to see if there is a match between KO# and JGI gene ID#
  for (j in 1:length(names(kegg_list)) ) {
    count_pathways_match_temp[j] = 0 # initialize count of number of pathways to which gene was assigned
    
    if (sum(kegg_list[[j]] == k_num) == 1){              # there is a match -- enter JGI ID into the new kegg pathway list
      kegg_list_jgi[[j]] = append(kegg_list_jgi[[j]], paste(jgi_num))
      count_pathways_match_temp[j] = count_pathways_match_temp[j] + 1
      #  }  else if (sum(kegg_list[[j]] == k_num) == 0) {          # the gene does not belong to any pathway
      #    no_pathway_match_gene = append(no_pathway_match_gene, gene_name)
    } else if (sum(kegg_list[[j]] == k_num) > 1) {       # there's something wrong if there is more than 1 match within a KEGG Pathway
      print(paste('Repeated KO# detected in KEGG PATHWAY: ', kegg_list[[j]][1], sep = ' '))
    }
  }
  
  if (sum(count_pathways_match_temp) == 0) { # gene had no match in any pathway
    no_pathway_match_gene = rbind(no_pathway_match_gene, c(gene_name, jgi_num, k_num))
  }
  
  # count the number of pathways to which gene was assigned
  count_pathways_match[z] = sum(count_pathways_match_temp)
  
}  

# housekeeping
no_ko_gene = no_ko_gene[-c(1),] # remove the first element, which is NA
no_pathway_match_gene = no_pathway_match_gene[-c(1),] # delete row 1 which is NA

#------- gene counts for methods section

# count the number of NA == # of genes with no KO assignment
sum(is.na(count_pathways_match))
length(no_ko_gene[,1]) # should be the same as above
sum(no_ko_gene[,1] == 'hypothetical protein')

# count the number of genes that did not match any KEGG pathways
length(no_pathway_match_gene[,1]) # length of 1 column
head(no_pathway_match_gene)
sum(no_pathway_match_gene[,1] == 'hypothetical protein') # should not be many

# save the genes with no K# assignments
write.csv(no_ko_gene,'dataset/jgi_downloaded/genes_with_no_KO_assignments.csv') 
# save the genes with no KEGG PATHWAY assignemnt
write.csv(no_pathway_match_gene,'dataset/jgi_downloaded/genes_with_no_KEGG_PATHWAY_assignments.csv') 

# distribution of genes being assigned
hist(count_pathways_match[count_pathways_match>0]) # only include non-zero assignments
max(count_pathways_match, na.rm = TRUE)
sum(count_pathways_match == 1,na.rm=TRUE) # majority of genes had only 1 KEGG assigments

# investigate the genes that were assigned many KEGG PATHWAYS
index = which(count_pathways_match == 0)
for (i in 1:length(index)){
  print(paste(data_jgi[[index[i]]][gene_name_index], ' (jgi ID #',data_jgi[[index[i]]][jgi_num_index],')',sep=''))
}

# count the number of KO number and JGI gene ID
# not really useful because there will always be a mismatch between number of JGI gene ID and KO numbers due to repeats
pathway_missing_genes = matrix() # pathways of missing genes
for (j in 1:length(names(kegg_list)) ) {    # loop through each KEGG pathway
  pathway = matrix() # initialize
  pathway = kegg_list_jgi[[j]]
  pathway = pathway[2:length(pathway)] # get rid of the first element, which is pathway description
  if (!sum(startsWith(pathway,'K')) == sum(startsWith(pathway,'6'))) {
    difference = sum(startsWith(pathway,'K')) - sum(startsWith(pathway,'6'))
    pathway_missing_genes[j] = print(paste(kegg_list_jgi[[j]][1], ' is missing ', difference , 'genes'))
  }
}
#------- END: gene counts for methods section




# .gmt file for GSEA
  # get rid of KO numbers from the KEGG Pathway gene list.
  # MAKE SURE that there is not already a file with the save file name. Otherwise it would append.
kegg_list_jgi_clean = kegg_list_jgi # with no KO numbers, only JGI gene IDs, JUST for counting below
for (z in 1:length(kegg_list_jgi)) {
  
  number = names(kegg_list_jgi)[z] # KEGG PATHWAY number only
  pathway = matrix()
  pathway = kegg_list_jgi[[z]] # KEGG PATHWAY and KO# and JGI gene ID#
  
  # indexes of JGI gene ID (and pathway description)
  index_to_keep = append(1, which(startsWith(pathway,'6'))) # always keep the first element (pathway description)
  
  # append KEGG pathway number to pathway descriptor
  path_descriptor = paste(number,kegg_list_jgi[[number]][1], sep = '_') 
  
  # SAVE to .gmt file
 # sink('gsea/gsea_input/vcor_kegg_path_jgi_gene_id.gmt',append=TRUE)
#    cat(paste(c(path_descriptor,pathway[index_to_keep],sep=" "),collapse='\t')) # append KEGG pathway number to pathway descriptor
#    cat('\n')                                                                   # create new line for next pathway
 # sink()
  
  # just for counting, clean up the list (only keep JGI gene IDs)
  kegg_list_jgi_clean[[z]] = kegg_list_jgi_clean[[z]][index_to_keep]
}
kegg_list_jgi = kegg_list_jgi_clean # for saving
save(kegg_list_jgi, file='dataset/kegg_downloaded/kegg_pathways_jgi.Rdata')

# print lengths for METHODS section 
num_k_assigned = lengths(kegg_list_jgi_clean,use.names=T) - 1 # number of K numbers in each pathway gene set (subtract the pathway name)
which.max(num_k_assigned) # identify the pathway gene set with the most number of K numbers assigned

# save the number of genes in each pathway (not perfect)
sink('dataset/kegg_downloaded/number_genes_per_kegg_pathway.txt')
  sort(num_k_assigned) # sort by number of K numbers assigned
sink()


