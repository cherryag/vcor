# Generate plots of time points vs. normalized expression of single genes

# Gene counts normalization method: normalize counts by size factor, i.e. library size (counts(dds, normalized=TRUE) )

#----- DESeq2 normalized gene counts
# dds <- estimateSizeFactors(dds) (perform on fresh dds)
# normalized_counts = counts(dds, normalized=TRUE)
  # This is just dividing each column of counts(dds) by sizeFactors(dds)
  # i.e. normalizing by library size
# It is possible that longer length genes recruit more reads, but because each gene is the same 
# across experimental conditions (t10m, t0c etc.), counts were NOT normalized by gene length.
# Keeping this in mind (i.e. not normalizing by gene length), I should not compare counts between genes (i.e.,
# "if a gene has length L and another with length 2L, you would also expect the second gene to have normalized
# counta that were twice as large." - Michael Love on a forum.)
#--------------------

# CREATED : 2/7/2020 by Cherry Gao
# UPDATED : 3/23/2020

# install.packages("readxl")

# libraries
library(DESeq2)
library(readxl)
library(ggplot2)
library(reshape2)


#---------- setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
save_dir='/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/virulence_genes/timelapse_plots' # where to save timelapse plots

load(file='deseq/dds_group_design.RData') # dds from DESeq2

load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data
load(file='dataset/virulence_genes_curated/jgi_gene_id_virulence_keywords.Rdata') # load virulence genes that are searched for by keywords
load(file='dataset/virulence_genes_curated/jgi_gene_id_virulence_manual_list.Rdata') # load virulence genes that were manually curated



# list of genes to plot (MANUALLY iterate)
#method_of_gene_curation = c('KO', 'KO_COG_prod','prod','manual')
#gene_list = paste0('jgi_id_and_name_',method_of_gene_curation[1])
gene_list = list()
# gene_list = jgi_id_and_name_KO
# gene_list = jgi_id_and_name_KO_COG_prod
# gene_list[['hemolysin']] = jgi_id_and_name_prod$hemolysin
# gene_list[['manual']] = jgi_id_and_name_manual # to make for loops go smoothly
# gene_list[['timelapse_gene']] = cbind("647169781", "cell division protein FtsZ")
gene_list[['timelapse_gene']] = cbind("647169941", "RNA polymerase, sigma 70 subunit, RpoD")



### extract normalized counts from DESeq2
# dds <- estimateSizeFactors(dds) # normalization only accounts for library size (not gene length) s.t. I can compare across samples
# sizeFactors(dds) # check size factor for each sample (noticed that samples with larger read counts lead to larger size factors, e.g. rep1TiC)
norm_c = as.data.frame(counts(dds, normalized=TRUE)) # "counts(dds, normalized=TRUE) provides counts scaled by size or normalization factors

## extract gene ID and gene product name from JGI database (not the same order as dds)
jgi_gene_id_column = 1
product_name_column = 5
jgi_gene_list = unlist(lapply(data_jgi,'[[',jgi_gene_id_column)) # extract gene ID numbers only (not in the same order as norm_c)
prod_list = unlist(lapply(data_jgi,'[[',product_name_column)) # extract gene product descriptions only

## FOR LOOP: make a line plot for each gene of interest
for (z in 1:length(gene_list)){ # for each category
  # make a new folder for each category
  # save_folder = paste0('virulence_genes/timelapse_plots/',names(gene_list)[z],'/sigma factors (Kimes 2012)')
  save_folder = paste0('timelapse_genes/plots')
  # dir.create(save_folder) # create the folder
    
  for (i in 1:length(gene_list[[z]][,1])){ # for each gene
#  for (i in 57:65){ # for each gene 
    
    rm(cont_data,muc_data,c1,c2,c3,m1,m2,m3) # initialize
    
    ## find the index of gene of data
    
    jgi_gene_id_of_interest = gene_list[[z]][i,1] # '647169116'
    jgi_gene_name_of_interest = gene_list[[z]][i,2] # 'OmpU'
    
    gene_idx = which(rownames(norm_c) == jgi_gene_id_of_interest) # index in dds
    
    ## for plotting, make small data frames for each gene of interest
    c1 = c(norm_c$rep1TIC[gene_idx], norm_c$rep1T10C[gene_idx], norm_c$rep1T60C[gene_idx])
    c2 = c(norm_c$rep2TIC[gene_idx], norm_c$rep2T10C[gene_idx], norm_c$rep2T60C[gene_idx]) 
    c3 = c(norm_c$rep3TIC[gene_idx], norm_c$rep3T10C[gene_idx], norm_c$rep3T60C[gene_idx])
    m1 = c(norm_c$rep1TIM[gene_idx], norm_c$rep1T10M[gene_idx], norm_c$rep1T60M[gene_idx])
    m2 = c(norm_c$rep2TIM[gene_idx], NaN,                       norm_c$rep2T60M[gene_idx]) # rep2T10M is missing sample
    m3 = c(norm_c$rep3TIM[gene_idx], norm_c$rep3T10M[gene_idx], norm_c$rep3T60M[gene_idx])
    time = c(0,10,60)
    rep = c(1,1,1,2,2,2,3,3,3)
    
    # prepare control data
    cont_data = cbind.data.frame(time,c1,c2,c3)
    cont_data = melt(cont_data,id='time')
      treatment = rep('control',length(rep))
    cont_data = cbind.data.frame(cont_data,rep,treatment)
  
    # prepare mucus data
    muc_data = cbind.data.frame(time,m1,m2,m3)
    muc_data = melt(muc_data,id='time')
      treatment = rep('mucus',length(rep))
    muc_data = cbind.data.frame(muc_data,rep,treatment)
    
    data = rbind(cont_data,muc_data)
    
  
    # CREATE dot plot (no lines)
    q = ggplot(data, aes(x=time, y=value, color=treatment)) + 
      geom_point(size=5, alpha=0.7, stroke=0, na.rm=TRUE) + 
      # (values=c('black', 'pink'))  # specify colors of the dots
      scale_color_manual(values = c('#969696', '#6baed6') ) 
    
    # labels and background color
    q + theme_bw() + labs(title=paste('expression over time of JGI gene ID',jgi_gene_id_of_interest),
                          subtitle=paste0('JGI database product name: ', prod_list[which(jgi_gene_list==jgi_gene_id_of_interest)],'\n', # name of gene
                                          paste('Literature search annotation name:', jgi_gene_name_of_interest) ),
                          x='time (minutes)', 
                          y='counts (normalized by library size)') 
    
  
    # SAVE figure (png first, for inspection)
    save_fn = paste0(save_folder,'/time_vs_counts_',jgi_gene_id_of_interest)
    ggsave(paste0(save_fn,'.pdf'), dpi=300, dev='pdf', height = 5.4, width = 6, units = "in") 
    ggsave(paste0(save_fn,'.png'), dpi=300, dev='png', height = 5.4, width = 6, units = "in") 
    
    
  
  }  # end of loop through each gene
} # end of for loop through each category of genes
