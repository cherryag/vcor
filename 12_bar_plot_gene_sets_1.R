# Bar plots of log2FC of set of genes
# Adapted from old script named log2plot.R from 2016

# CREATED : 1/30/2020 by Cherry Gao
# UPDATED : 2/9/2020  added virulence gene set (manual)

# 1) specify variable 'i' which is the indexes in data_concat of genes in the gene set that is desired
# 2) create bar plot of foldchange for each gene set



############### SET UP ###############
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/gsea')
load(file='duplicate_knum_histograms/duplicate_knum_all.Rdata') # K numbers with multiple JGI gene ID assigments
load(file='kegg_downloaded/kegg_pathways.Rdata') # KEGG pathway - K number assignment
load(file='deseq_kegg_concat/60muc_vs_60cont_alpha_0.05_KEGGconcat.Rdata') # load data_concat
  # housekeeping on data_concat
  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
  data_concat$padj = as.numeric(paste(data_concat$padj))

  # manually change pairwise comparisons
  pairwise_title = '60m vs 60c' # for figure titles
  pairwise_fn = '60m_60c_' # for save file names
  
# initialize
geneset_i = list()

# libraries
library(ggplot2)
library(stringr)

############### PROCESS GENE SETS ###############

##### VIRULENCE GENES by key words (added 2/9/2020 manually)
description = 'Virulence genes by key words'
virulence_keywords = c('','','') # MANUAL LIST - keywords for virulence genes (in excel sheet, 'virulence_genes.xls')

# find virulence genes
virulence_i = vector() # initialize
for (z in 1:length(virulence_keywords)){
  keyword = virulence_keywords[z]
  virulence_i = c(grep(keyword, data_concat$Product_name,ignore.case=T), 
                 grep(keyword, data_concat$geneProduct,ignore.case=T)) 
}

# cleanup identified virulence genes
virulence_i = unique(virulence_i) # get rid of duplicate indexes
data_concat$jgi_gene_id[virulence_i] # check JGI gene ID ordering
geneset_i[[description]] = virulence_i # FINALIZE

##### gene set 1: k numbers with multiple assignments of JGI gene IDs
for (rep_n in c(5,6,7,20,50)) { # number of JGI gene IDs assigned to k number
  k_i = which(knum_to_check$n==rep_n)
  k_num = knum_to_check$KO_new[k_i] # print k numbers
  k_num_description = 
  
  for (j in 1:length(k_num)){ # for each k num with a certain number of JGI gene IDs assigned
    i_temp = grep(k_num[j], str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*")) # extract index in data_concat
    description = paste(k_num[j],'(', data_concat$geneProduct[i_temp[1]],'):', rep_n, 'JGI gene ID assignments')
    geneset_i[[description]] = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
  }
}


##### gene set 2: RIBOSOME genes by keywords 
description = 'RIBOSOME genes by keywords "rRNA","ribosom","[J] (COG category)"'
# find
ribosome_i = c(grep('rRNA',data_concat$Product_name), 
               grep('rRNA',data_concat$geneProduct),
               grep('ribosom',data_concat$Product_name),
               grep('ribosom',data_concat$geneProduct),
               grep('[J]',data_concat$COG_category)) # J category = 'Translation, ribosomal structure and biogenesis'
# order by JGI gene ID
ribosome_i = unique(ribosome_i) # get rid of duplicate indexes
  ribosome_i = ribosome_i[order(data_concat$jgi_gene_id[ribosome_i])] # order i by JGI gene ID (to identify operons)
  data_concat$jgi_gene_id[ribosome_i] # check JGI gene ID ordering
  plot(1:length(data_concat$jgi_gene_id[ribosome_i]),data_concat$jgi_gene_id[ribosome_i]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = ribosome_i # FINALIZE


##### gene set 3: RIBOSOME genes by KEGG pathway categorization (top GSEA pathways)
description = 'RIBOSOME genes by KEGG (03011 + 03010)'
k_num_interest = vector() # initialize
# find
kegg_path_name_1 = "03011_Ribosome [BR:vct03011]" # from GSEA
  k_num_interest = append(k_num_interest, kegg_pathways[[kegg_path_name_1]])
  kegg_pathways[[kegg_path_name_1]] # check by printing
kegg_path_name_2 = "03010_Ribosome [PATH:vct03010]" # from GSEA
  kegg_pathways[[kegg_path_name_2]] # check by printing
  k_num_interest = append(k_num_interest, kegg_pathways[[kegg_path_name_2]])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
  data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
  plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE

# investigate intersection of the two ribosome gene sets
geneset_i[['RIBOSOME genes by KEGG (03011 + 03010)']] %in% geneset_i[['RIBOSOME genes by keywords "rRNA","ribosom","[J] (COG category)"']]


##### gene set 4: RIBOSOME BIOGENESIS genes by KEGG
description = '03009_Ribosome biogenesis [BR:vct03009]' # FROM GSEA
k_num_interest = vector() # initialize
# find
kegg_pathways[[description]] # make sure the name exists
k_num_interest = append(k_num_interest, kegg_pathways[[description]])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 5: amino acid METABOLISM genes by KEGG (no ordering)
description = 'amino acid metabolism genes by KEGG (00230 + 00340 + 00920 + 00260 + 00620)'
k_num_interest = vector() # initialize
# find
k_num_interest = append(k_num_interest, kegg_pathways[['00230_Purine metabolism [PATH:vct00230]']])
k_num_interest = append(k_num_interest, kegg_pathways[['00340_Histidine metabolism [PATH:vct00340]']])
k_num_interest = append(k_num_interest, kegg_pathways[['00260_Glycine, serine and threonine metabolism [PATH:vct00260]']])
k_num_interest = append(k_num_interest, kegg_pathways[['00620_Pyruvate metabolism [PATH:vct00620]']])
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 6: SULFUR METABOLISM genes by KEGG
description = '00920_Sulfur metabolism [PATH:vct00920]' # FROM GSEA
k_num_interest = vector() # initialize
# find
kegg_pathways[[description]] # make sure the name exists
k_num_interest = append(k_num_interest, kegg_pathways[[description]])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
  data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
  plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 7: late-stage protein synthesis genes by KEGG (no ordering)
description = 'late-stage protein synthesis genes by KEGG (03110 + 03029 + 03016)'
k_num_interest = vector() # initialize
# find
k_num_interest = append(k_num_interest, kegg_pathways[['03110_Chaperones and folding catalysts [BR:vct03110]']])
k_num_interest = append(k_num_interest, kegg_pathways[['03029_Mitochondrial biogenesis [BR:vct03029]']])
k_num_interest = append(k_num_interest, kegg_pathways[['03016_Transfer RNA biogenesis [BR:vct03016]']])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 8: PHOTOSYNTHESIS genes (??) by KEGG
description = '00194_Photosynthesis proteins [BR:vct00194]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 9: MOTILITY genes by KEGG
description = '02035_Bacterial motility proteins [BR:vct02035]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 10: FLAGELLA genes by KEGG
description = '02040_Flagellar assembly [PATH:vct02040]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 11: CHEMOTAXIS genes by KEGG
description = '02030_Bacterial chemotaxis [PATH:vct02030]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 12: SECRETION genes by KEGG (no ordering)
description = 'SECRETION genes by KEGG (03070 + 02044 + 03060)'
k_num_interest = vector() # initialize
# find
k_num_interest = append(k_num_interest, kegg_pathways[['03070_Bacterial secretion system [PATH:vct03070]']])
k_num_interest = append(k_num_interest, kegg_pathways[['02044_Secretion system [BR:vct02044]']])
k_num_interest = append(k_num_interest, kegg_pathways[['03060_Protein export [PATH:vct03060]']])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 13: SIGNAL TRANSDUCTION genes by KEGG (no ordering)
description = 'SIGNALING PROTEINS + KINASES genes by KEGG (99995 + 01001)'
k_num_interest = vector() # initialize
# find
k_num_interest = append(k_num_interest, kegg_pathways[['99995_Signaling proteins']])
k_num_interest = append(k_num_interest, kegg_pathways[['01001_Protein kinases [BR:vct01001]']])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 14: BIOFILM genes by KEGG
description = '05111_Biofilm formation - Vibrio cholerae [PATH:vct05111]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 15: VIRULENCE genes by KEGG (no ordering)
description = 'VIRULENCE genes by KEGG (01501 + 01002 + 02042 + 02022 + 02020)'
k_num_interest = vector() # initialize
# find
k_num_interest = append(k_num_interest, kegg_pathways[['01501_beta-Lactam resistance [PATH:vct01501]']])
k_num_interest = append(k_num_interest, kegg_pathways[['01002_Peptidases and inhibitors [BR:vct01002]']])
k_num_interest = append(k_num_interest, kegg_pathways[['02042_Bacterial toxins [BR:vct02042]']])
k_num_interest = append(k_num_interest, kegg_pathways[['02022_Two-component system [BR:vct02022]']])
k_num_interest = append(k_num_interest, kegg_pathways[['02020_Two-component system [PATH:vct02020]']])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### gene set 16: QUORUM SENSING genes by KEGG
description = '02024_Quorum sensing [PATH:vct02024]'
k_num_interest = vector() # initialize
# find
k_num_interest = kegg_pathways[[description]]
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE

# SAVE
save(geneset_i,file=paste0('gene_set_bar_graph/genesets_',pairwise_fn,'for_bar.Rdata'))




##### MODIFIED gene sets (RE-setup)
# rm(list = ls()) # clear environment
# setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/gsea')
# load(file='duplicate_knum_histograms/duplicate_knum_all.Rdata') # K numbers with multiple JGI gene ID assigments
# load(file='kegg_downloaded/kegg_pathways.Rdata') # KEGG pathway - K number assignment
# load(file='deseq_kegg_concat/10muc_vs_10cont_alpha_0.05_KEGGconcat.Rdata') # load data_concat
  # housekeeping on data_concat
#  data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
#  data_concat$padj = as.numeric(paste(data_concat$padj))

# load(file='gene_set_bar_graph/genesets_for_bar.Rdata')

# define 'not in' function
'%!in%' <- function(x,y)!('%in%'(x,y))

##### rRNA keyword only
description = 'rRNA only'
# find
i_temp = c(grep('rRNA',data_concat$Product_name), 
           grep('rRNA',data_concat$geneProduct))
# order by JGI gene ID
i_temp = unique(i_temp) # get rid of duplicate indexes
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
  data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE


##### bacterial CHEMOTAXIS genes WITHOUT MCPs
mcp_i = geneset_i$`K03406 ( methyl-accepting chemotaxis protein ): 50 JGI gene ID assignments`
chemotaxis_i = geneset_i$`02030_Bacterial chemotaxis [PATH:vct02030]`
chemotaxis_i_new = chemotaxis_i[chemotaxis_i %!in% mcp_i]
description = '02030_Bacterial chemotaxis genes WITHOUT MCPs'
  i_temp = unique(chemotaxis_i_new) # get rid of duplicate indexes
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE

##### bacterial MOTILITY genes WITHOUT MCPs
mcp_i = geneset_i$`K03406 ( methyl-accepting chemotaxis protein ): 50 JGI gene ID assignments`
motility_i = geneset_i$`02035_Bacterial motility proteins [BR:vct02035]`
motility_i_new = motility_i[motility_i %!in% mcp_i]
description = '02035_Bacterial motility proteins WITHOUT MCPs'
  i_temp = unique(motility_i_new) # get rid of duplicate indexes
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE

##### bacterial MOTILITY genes WITHOUT MCPs and WIHTOUT FLAGELLA genes
mcp_i = geneset_i$`K03406 ( methyl-accepting chemotaxis protein ): 50 JGI gene ID assignments`
flagella_i = geneset_i$`02040_Flagellar assembly [PATH:vct02040]`
motility_i = geneset_i$`02035_Bacterial motility proteins [BR:vct02035]`
chemotaxis_i = geneset_i$`02030_Bacterial chemotaxis [PATH:vct02030]`

motility_i_new = motility_i[motility_i %!in% c(mcp_i, flagella_i, chemotaxis_i)]
description = '02035_Bacterial motility proteins WITHOUT MCPs, WIHTOUT FLAGELLA genes, WITHOUT CHEMOTAXIS genes'
  i_temp = unique(motility_i_new) # get rid of duplicate indexes
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE



##### Two-component system
description = 'TWO-COMPONENT system by KEGG (02022 + 02020)'
k_num_interest = vector() # initialize
# find
kegg_path_name_1 = '02022_Two-component system [BR:vct02022]' # from GSEA
  k_num_interest = append(k_num_interest, kegg_pathways[[kegg_path_name_1]])
  kegg_pathways[[kegg_path_name_1]] # check by printing
kegg_path_name_2 = "02020_Two-component system [PATH:vct02020]" # from GSEA
  kegg_pathways[[kegg_path_name_2]] # check by printing
  k_num_interest = append(k_num_interest, kegg_pathways[[kegg_path_name_2]])
# order by JGI gene ID
k_num_interest = unique(k_num_interest) # get rid of duplicates
  i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
  i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])] # order i by JGI gene ID (to identify operons)
data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
geneset_i[[description]] = i_temp # FINALIZE

#####  for-loop for single KEGG categories
description_list = c('00230_Purine metabolism [PATH:vct00230]',    # amino acid metabolism separately
                     '00340_Histidine metabolism [PATH:vct00340]',
                     '00260_Glycine, serine and threonine metabolism [PATH:vct00260]',
                     '00620_Pyruvate metabolism [PATH:vct00620]',
                     '03110_Chaperones and folding catalysts [BR:vct03110]', # late-stage protein synthesis genes by KEGG (no ordering)
                     '03029_Mitochondrial biogenesis [BR:vct03029]',
                     '03016_Transfer RNA biogenesis [BR:vct03016]',
                     '03070_Bacterial secretion system [PATH:vct03070]', # secretion
                     '02044_Secretion system [BR:vct02044]',
                     '03060_Protein export [PATH:vct03060]',
                     '01501_beta-Lactam resistance [PATH:vct01501]', # virulence
                     '01002_Peptidases and inhibitors [BR:vct01002]',
                     '02042_Bacterial toxins [BR:vct02042]')

for (j in 1:length(description_list)) {
  description = description_list[j]
  k_num_interest = vector() # initialize
  # find
  k_num_interest = append(k_num_interest, kegg_pathways[[description]])
    i_temp = which(str_extract(data_concat$KO, "(?<=KO:)[0-9A-Za-z]*") %in% k_num_interest) # some of the k numbers assigned to pathways (from another Vcor strain than BAA) do not exist in my Vcor's genome
    i_temp = i_temp[order(data_concat$jgi_gene_id[i_temp])]
  data_concat$jgi_gene_id[i_temp] # check JGI gene ID ordering
  plot(1:length(data_concat$jgi_gene_id[i_temp]),data_concat$jgi_gene_id[i_temp]) # visualize arrangement (by JGI gene ID)
  geneset_i[[description]] = i_temp # FINALIZE

}


# RE-SAVE
save(geneset_i,file=paste0('gene_set_bar_graph/genesets_',pairwise_fn,'for_bar_modified.Rdata'))

############### END OF PROCESS GENE SETS ###############




############### BAR PLOT ###############
# rm(list = ls()) # clear environment
# setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/gsea')
# load(file='duplicate_knum_histograms/duplicate_knum_all.Rdata') # K numbers with multiple JGI gene ID assigments
# load(file='kegg_downloaded/kegg_pathways.Rdata') # KEGG pathway - K number assignment
# load(file='deseq_kegg_concat/60muc_vs_60cont_alpha_0.05_KEGGconcat.Rdata') # load data_concat
# load(file='gene_set_bar_graph/genesets_60m_60c_for_bar_modified.Rdata') 
# housekeeping on data_concat
# data_concat$foldchange = as.numeric(paste(data_concat$foldchange))
# data_concat$padj = as.numeric(paste(data_concat$padj))
# manually change pairwise comparisons
# pairwise_title = '60m vs 60c' # for figure titles
# pairwise_fn = '60m_60c_' # for save file names

# threshold on padj for different colors of bars according to significance
p_thresh = 0.05 

# loop through all repeated k numbers
# for (z in 1:length(geneset_i)) {
# for (z in z_trouble) {
for (z in 1:length(geneset_i)) {
  title_text = paste(pairwise_title,'\n',names(geneset_i)[z],'\n padj threshold = ',p_thresh,'\n grey = n.s.; red = significant')
  
  i = geneset_i[[names(geneset_i)[z]]] # extract indexes in data_concat of genes to be plotted as bars
  data_subset = data_concat[i,]
  data_subset$bar_labels = paste0(data_concat$geneProduct[i],"_", data_concat$jgi_gene_id[i])
  
  # ggplot : bar graph of all genes in the set
  g = ggplot(data_subset)  # needs data frame
  g = g + geom_bar( aes(x=bar_labels, y=foldchange, fill=(padj < p_thresh) ), stat="identity") # specify data + labels to plot
  g = g + coord_flip()  # horizontal bars
  g = g + scale_fill_manual(values = c('grey', 'red') )  # false, true bar colors
  g = g + ggtitle(title_text) + ylab('log2 fold change') + xlab('gene product _ JGI gene ID')
  g = g + theme(plot.title = element_text(hjust = 1, size = 10, color='black'))

  print(g)
  
  # save figure
  if (z == 16) { # trouble for some reason
    text = str_replace(names(geneset_i)[z],'/signal transduction systems component','') # for z = 16, K# with 20 JGI gene ID assigments
    save_fn = paste0('geneset_bars_fc_',pairwise_fn,text,'.eps')
  }
  else {
    save_fn = paste0('geneset_bars_fc_',pairwise_fn,names(geneset_i)[z],'.eps')
  }
  ggsave(save_fn, dpi=300, dev='eps') # NOTE TO SELF: make the font smaller so sizes in inches make sense
  # ggsave(save_fn, dpi=300, dev='eps', height = 10, width = 15, units = "in") # NOTE TO SELF: make the font smaller so sizes in inches make sense
}


# after inspecting plots, manually picked out figures that had trouble and should be re-done
# z_trouble = c(16,18,21,23,24,25,27,28,31) 

        
############### END OF BAR PLOT ###############


        



        
###### OLD

# other genes that may be of interest, from OLD code (2016)
# ABC transporter, t10m vs. t10c
# cGMP (rows 56:60)
# che genes (rows 61:71)
# chitin (rows 72:107)
# flagella (rows 107:185)
# phosphate (rows 186:368)
# pilus (rows 369:414)
# quorum sensing (rows 415:427)
# ribosome (rows=428:515)
# sodium (rows 516:553)
# virulence (rows 554:560)


### R code (works) : bar graph of all genes in the set
# par(mar=c(5, 12, 4, 3)) # bottom, left, top and right margin sizes (replace mar with mai for inches)
# par(mai=c(1,1,1,1), mgp=c(3,1,0)) # mgp = axis title, labels, and line margin lines (but applies to both x and y axes)
# barplot(data_concat$foldchange[i], 
#         col = ifelse(data_concat$padj[i]<p_thresh, "red","grey"),
#         horiz = TRUE, 
#         names.arg = as.character(paste0(data_concat$geneProduct[i],"_", data_concat$jgi_gene_id[i])),
#         font.axis=1, las=2, # label directions
#         cex.names=0.5, cex.lab=1, cex.axis=1, cex.main=1, # label font sizes
#         xlab = "log2 fold change", main = title_text)
  