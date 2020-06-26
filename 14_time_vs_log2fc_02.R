# Generate plots of time points vs. log2FC of single genes (i.e. virulence genes)

# CREATED : 3/25/2020 by Cherry Gao
# UPDATED : 4/18/2020  (to _02.R)
#           Plot quorum sensing genes (from Kimes et al. 2012) -- make final figures
#           Need to improve to become more versatile and less manual (right now, need to manually iterate for different gene groups). 


# libraries
library(ggplot2)
library(reshape2)
library(tidyverse)


#---------- setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
save_dir='/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/virulence_genes/timelapse_log2fc_plots' # where to save timelapse plots
  # dir.create(save_dir)

load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data
# load(file='dataset/virulence_genes_curated/jgi_gene_id_virulence_manual_list_trimmed.Rdata') # load virulence genes that were manually curated
# load(file='dataset/virulence_genes_curated/jgi_gene_id_virulence_manual_list_trimmed_ncbi.Rdata') # load virulence genes that were manually curated after cross-check with NCBI [ADDED 4/25/2020]
# load(file='dataset/virulence_genes_curated/jgi_gene_id_sets_kimes_qs_v2.Rdata') # load quorum sensing genes (Kimes et al.)
# load(file='dataset/virulence_genes_curated/jgi_gene_id_sets_kimes_protease.Rdata') # load protease genes (Kimes et al.)
load(file='dataset/virulence_genes_curated/jgi_gene_id_vps.Rdata') # VPS biofilm genes

# list of genes to plot (MANUALLY iterate)
# method_of_gene_curation = c('KO', 'KO_COG_prod','prod','manual')
# gene_list = paste0('jgi_id_and_name_',method_of_gene_curation[1])
gene_list = list()
# gene_list = jgi_id_and_name_KO
# gene_list = jgi_id_and_name_KO_COG_prod
# gene_list[['hemolysin']] = jgi_id_and_name_prod$hemolysin
# gene_list[['virulence genes']] = jgi_id_and_name_manual # to make for loops go smoothly
# gene_list = jgi_gene_id_set_kimes
# gene_list[['ncbi_checked']] = jgi_id_and_name_manual # the [[1]] makes it easier for later
gene_list[['vps biofilm genes']] = jgi_id_and_name_manual # the [[1]] makes it easier for later

#----- for QS genes-------
# luxM = cbind('647173492','N-(3-hydroxybutanoyl)-L- homoserine lactone synthase LuxM') # not in Kimes
# gene_list[['quorum sensing kimes']] = rbind(jgi_gene_id_set_kimes$`Quorum sensing (added)`,luxM)
#  geneset_description = 'quorum sensing (Kimes et al)'
#-------------------------


# sequentially load and store log2FC 10m vs. 10c, and 60m vs. 60c DESeq2 results
# 10 min
load(file='deseq/deseq_kegg_concat/10muc_vs_10cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_10 = data_concat
rm(data_concat)
# 60 min
load(file='deseq/deseq_kegg_concat/60muc_vs_60cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_60 = data_concat
rm(data_concat)

## extract gene ID and gene product name from JGI database (not the same order as dds)
# jgi_gene_id_column = 1
# product_name_column = 5
# jgi_gene_list = unlist(lapply(data_jgi,'[[',jgi_gene_id_column)) # extract gene ID numbers only (not in the same order as norm_c)
# prod_list = unlist(lapply(data_jgi,'[[',product_name_column)) # extract gene product descriptions only

# [OPTIONAL] MANUALLY determine indexes of genes to plot on line graph
gene_list # visualize 
#----- for VIRULENCE genes (paired only with 'jgi_id_and_name_manual'):
# gene_index_to_plot_list = list(setdiff(c(1:5,7:9,12,15:20),c(2,5,7,12,15:20)),  c(1:5,7:9,12,15:20),c(2,5,7,12,15:20),c(6,13,14), c(1:20), c(21:27),c(21,24,25),
#                                       c(2,8,16,17,18,20), c(1,5,7,9,12,15)) # manual, indexes of genes
# geneset_description_list = list('Vcor virulence genes (no toxins)',
#                                 'Vcor virulence genes',
#                                'toxins',
#                                 'antibiotic resistance',
#                                 'all virulence genes',
#                                 'all zinc metalloproteases',
#                                 'select zinc metalloproteases',
#                                 'positively regulated virulence genes',
#                                 'negatively or not signficantly regulated virulence genes') # manual, for plot title
#-------------------------

#----- for VIRULENCE genes (paired only with 'jgi_gene_id_virulence_manual_list_trimmed_ncbi.Rdata'):
# gene_index_to_plot_list = list(c(1:5), c(6:10), c(11:14), c(15:21), c(22:23), c(24:31), c(1,32:35), c(1,35), c(37:52)) # manual, indexes of genes
# geneset_description_list = list('transcriptional regulators',
#                                 'host colonization',
#                                 'biofilm matrix',
#                                 'host damage + toxins',
#                                 'antibiotic resistance',
#                                 'type 6 secretion system',
#                                 'quorum sensing inducers',
#                                 'qs master regulators',
#                                 'putative vps region') # manual, for plot title
#-------------------------


#----- for QS genes:------
# gene_index_to_plot_list = list(c(1), c(2,3,4), c(11,12),c(5,6,7,8,9,10),c(1,4,13))
# geneset_description_list = list('CAI-1','AI-2','AI-3','general QS','CAI-1, AI-1, AI-2 producers')
#-------------------------

#----- when there's only 1 category
gene_index_to_plot_list = list(c(1:length(gene_list[[1]][,1])))
geneset_description_list = list(names(gene_list))
#------------------------

for (geneset_to_plot in c(1)){  # previously, manually and iteratively change
  
gene_index_to_plot = gene_index_to_plot_list[[geneset_to_plot]]
geneset_description = geneset_description_list[[geneset_to_plot]]

## FOR LOOP: make a line plot for each gene of interest
for (z in 1:length(geneset_description)){ # for each category
  
  # make a new folder for each category
  save_folder = paste0(save_dir,'/',names(gene_list)[z])
  # dir.create(save_folder) # create the folder
    
  #----- prepare dataset to plot
  
  # extract indexes of genes in data_concat 
  gene_idx_10 = match(as.numeric(unlist(gene_list[[1]][,1])), data_concat_10$seq_data_gene_id)
      gene_idx_10 = gene_idx_10[gene_index_to_plot]
  gene_idx_60 = match(as.numeric(unlist(gene_list[[1]][,1])), data_concat_60$seq_data_gene_id)
      gene_idx_60 = gene_idx_60[gene_index_to_plot]
  gene_name_list = unlist(gene_list[[1]][gene_index_to_plot,2])
  gene_name_list = paste0(unlist(gene_list[[1]][gene_index_to_plot,1]),'_',gene_name_list) # paste the JGI gene ID numbers
  
  # concatenate gene data to plot
  fc_10 = matrix(); padj_10 = matrix(); fc_60 = matrix(); padj_60 = matrix(); data_plot_temp = matrix() # initialize
    fc_10 = as.numeric(paste(data_concat_10$foldchange[gene_idx_10]))
    padj_10 = as.numeric(paste(data_concat_10$padj[gene_idx_10]))
  data_10 = cbind(10, fc_10, padj_10, as.character(gene_name_list))
  
    fc_60 = as.numeric(paste(data_concat_60$foldchange[gene_idx_60]))
    padj_60 = as.numeric(paste(data_concat_60$padj[gene_idx_60]))
  data_60 = cbind(60, fc_60, padj_60, as.character(gene_name_list))
    
  # FUTURE: should be the actual addition time (negative number)
  data_0 = cbind(0, 0, 0, as.character(gene_name_list)) # 0 time point -- assumed to have log2fc of 0
    
  # format gene data for plotting
  data_plot_temp = rbind(data_0,data_10,data_60)
  data_plot <- data.frame(
                      time = as.numeric(data_plot_temp[,1]),
                      log2fc = as.numeric(data_plot_temp[,2]),
                      padj = as.numeric(data_plot_temp[,3]),
                      gene = data_plot_temp[,4]
                       )
    
   #----- VERSION 1: line + dot plot
   g = ggplot(data_plot, aes(x=time, y=(log2fc))) + # for absolute fc, 2^log2fc
     #  geom_line(aes(color=gene)) + 
       geom_point(size=3,aes(color=gene, alpha=(padj<0.05))) + 
       geom_text(data=data_plot %>% filter(time == 10), aes(label = gene, # time==last(time)
                                                                  x = time + 15,  # time - 15
                                                                  y = log2fc + 0.05, 
                                                                  color = gene), size=2) # label each line 
#   g
    
    # labels and background color
    g = g + theme_bw() + 
        labs(title=paste('Log2 fold change of ',geneset_description),
                          subtitle='mucus vs. control comparison at each time point; t0 = time point before addition, log2fc assumed to be 1',
                          x='Time (minutes)', 
                          y='Log2 fold change') +
        theme(text=element_text(size=5)) + 
        theme(legend.position='top') #+ 
   #     ylim(-4,6)
   # ylim(-2.5,2.5)
       
    
    #---- VERSION 2: dot plot (used for protease, 4/20/2020)
    # data_plot_v2 = rbind(data_plot[data_plot$time == 10,],data_plot[data_plot$time == 60,]) # only 10 and 60 min data points
    # g = ggplot(data_plot_v2, aes(x=time, y=(2^log2fc))) + # for absolute fc, 2^log2fc
        # geom_line(aes(color=gene)) + 
    #     geom_point(size=3,aes(color=gene, alpha=(padj<0.01))) + 
    #     geom_text(data=data_plot %>% filter(time == last(time)), aes(label = gene, 
    #                                                               x = time - 15, 
    #                                                               y = 2^log2fc + 0.05, 
    #                                                               color = gene), size=2) # label each line 
    
    # line at fold-change of 1
    # g = g + geom_hline(yintercept=1, linetype="dashed", color = "black")
    
    # labels and background color
    # g = g + theme_bw() + 
    #   labs(title=paste('Log2 fold change of ',geneset_description),
    #       subtitle='mucus vs. control comparison at each time point; t0 = time point before addition, log2fc assumed to be 1',
    #       x='Time (minutes)', 
    #       y='absolute fold change') +
    #  theme(text=element_text(size=5)) + 
    #  theme(legend.position='top')
    
  
    #----- SAVE figure (png first, for inspection)
    save_fn = paste0(save_folder,'/time_vs_log2fc_',geneset_description)
    ggsave(paste0(save_fn,'_dots_v2.pdf'), dpi=300, dev='pdf', height = 6, width = 5.5, units = "in",useDingbats=FALSE) # useDingbats for ggplot -> Adobe Illustrator transition (markers are not fonts)
    ggsave(paste0(save_fn,'_dots_v2.png'), dpi=300, dev='png', height = 4, width = 5.5, units = "in") 
    
} # end of for loop through each category of genes

} 