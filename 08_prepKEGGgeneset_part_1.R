# Create geneset database defined by KEGG PATHWAY 
# PART 1 ONLY (clean up KEGG-downloaded database)

# Categorize genes (= KEGG Orthology (KO) numbers) for Vcor OCN014.
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


#----------PART 1: clean KEGG dataset----------
# Goal: create a spreadsheet with KEGG PATHWAY followed by columns of KO#s
# e.g.: KEGG PATHWAY; KO#; KO#; KO#...

# set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')

# read KEGG table *downloaded from internet)
file_name = 'dataset/kegg_downloaded/vct00001.txt'
raw_kegg = readLines(file_name) # character object; refer to each line as line = raw_kegg[line_number]
head(raw_kegg) # look at table

# BRITE hierarchies start at line 3012 with text, "A09180 Brite Hierarchies"
line_brite = 3012 # discerned from Excel sheet
raw_kegg[line_brite]

# there are categories that are "not included in Pathway or Brite"
  # gene sets like "unclassified metabolism" or "function unknown"
line_other = 5128 # from Excel
raw_kegg[line_other]

# range of lines to parse through (KEGG PATHWAY only)
line_range = 1:3011

# parse KEGG table
  # lines with letter "C" preceeds a KEGG PATHWAY gene set names (e.g. '00010 Glycolysis / Gluconeogenesis [PATH:vct00010]')
  # lines with letter "D" preceeds a gene (JV59_#####) and KO number in the next column (e.g. 'D      JV59_00740 rod shape-determining protein MreB	// K03569 mreB; rod shape-determining protein MreB and related proteins')
kegg_list = list() # initialize container

for (l in line_range){
  line = raw_kegg[l]
  # 1) extract KEGG PATHWAY category
  if(substr(line, 1, 1) == 'C') {                                          # if line begins with 'C'; function "substr(character_vector, start, stop)"
    a = unlist(strsplit(line, ' +'))                                        # function "strsplit(character_vector, "regular expressions to use for splitting")"
    path_number = a[[2]]                                                    # choose second strsplit element, containing BRITE category number
    path_name = gsub('\t', '', paste(a[3:length(a)], collapse=' '))         # extract name of BRITE category; function "gsub(pattern, replacement, vector_string)"
    # Append parsed KEGG data in container 'kegg_list'
    if (path_number %in% names(kegg_list)){stop('Error: duplicate names')}  # error if duplicate header category (starting with 'C')
    kegg_list[[path_number]] = c(path_name)
  
  # 2) extract KO number
  } else if(substr(line, 1, 1) == 'D') {
    a = unlist(strsplit(line, '\t'))[[2]]                           # parse by tab (i.e. column); then select 2nd tab field
    knumber = gsub(' .*', '', a)                                    # remove everything following space, leaving only number (in quotation marks)
#     knumber = gsub('\"', '', knumber)                             # remove quotation marks from the number, leaving only the number
    kegg_list[[path_number]] = c(kegg_list[[path_number]], knumber) # append KO number to kegg_list under the most recent path_number
  } else {
    next
  }
}

# Get rid of repeated KO numbers within each pathway [NEW as of 2/6/2020]
# There repeated KO numbers within each PATHWAY because several "JV####" genes can be assigned to one KO#. 
for (z in 1:length(kegg_list)){
  k_list = kegg_list[z]
  kegg_list[[z]] = unique(k_list[[1]])
}

# Extract and delete PATHWAY categories ("C" lines) with no assigned KO numbers ("D" lines)
kegg_list_empty = list()
for(number in names(kegg_list)){                     # names = BRITE category numbers
  if(length(kegg_list[[number]]) == 1){              # length == 1 means only contains the name of the category number
    kegg_list_empty[[number]] = kegg_list[[number]]  # extract and keep track of empty categories
    kegg_list[[number]] = NULL                       # delete empty categories from final KEGG list
  }
}

# Append category number to category name + KO numbers, then SAVE to file
sink('kegg_downloaded/vcor_kegg_pathways.txt')
for(number in names(kegg_list)){
  out = c(number, kegg_list[[number]])    # append category number to category name + KO numbers
  cat(paste(unique(out), collapse='\t'))  # print out only unique KO numbers, as tab-separated
  cat('\n')                               # create new line
}
sink()

# SAVE empty BRITE (or PATHWAY) categories
sink('kegg_downloaded/vcor_kegg_pathways_empty.txt')
for(number in names(kegg_list_empty)) {
  out = c(number, kegg_list_empty[[number]])    # append category number to category name + KO numbers
  cat(paste(unique(out), collapse='\t'))        # print out only unique KO numbers, as tab-separated
  cat('\n')                                     # create new line
}
sink()

# print lengths
length(kegg_list)
length(kegg_list_empty)
length(kegg_list) + length(kegg_list_empty)
num_k_assigned = lengths(kegg_list,use.names=T) # number of K numbers in each pathway gene set
which.max(num_k_assigned) # identify the pathway gene set with the most number of K numbers assigned
sort(num_k_assigned) # sort by number of K numbers assigned

# append pathway descriptors to KEGG pathway numbers (for save to .Rdata)
kegg_pathways = list()
for (j in 1:length(names(kegg_list))) {
  pathway_group = kegg_list[[j]]
  new_name = paste(names(kegg_list[j]), pathway_group[1], sep='_') # paste together KEGG pathway number and descriptor
  kegg_pathways[[new_name]] = kegg_list[[j]][2:length(kegg_list[[j]])] # only the KO numbers (element 1 = descriptor of pathway)
}
# save as .RData (so that I don't have to read tables again)
save.image(file='kegg_downloaded/kegg_pathways.Rdata')

#----------END of PART 1: clean KEGG dataset----------

