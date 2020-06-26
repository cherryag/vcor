# Generate plots of time points vs. log2FC of single genes (i.e. virulence genes)

# CREATED : 3/25/2020 by Cherry Gao


# libraries
library(ggplot2)
library(reshape2)


#---------- setup
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
save_dir='/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/virulence_genes/timelapse_log2fc_plots' # where to save timelapse plots
  # dir.create(save_dir)

load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data
load(file='dataset/virulence_genes_curated/jgi_gene_id_virulence_manual_list_trimmed.Rdata') # load virulence genes that were manually curated


# list of genes to plot (MANUALLY iterate)
# method_of_gene_curation = c('KO', 'KO_COG_prod','prod','manual')
# gene_list = paste0('jgi_id_and_name_',method_of_gene_curation[1])
gene_list = list()
# gene_list = jgi_id_and_name_KO
# gene_list = jgi_id_and_name_KO_COG_prod
# gene_list[['hemolysin']] = jgi_id_and_name_prod$hemolysin
gene_list[['virulence genes']] = jgi_id_and_name_manual # to make for loops go smoothly


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
#  prod_list = unlist(lapply(data_jgi,'[[',product_name_column)) # extract gene product descriptions only

# MANUALLY determine indexes of genes to plot on line graph
gene_list # visualize 
gene_index_to_plot_list = list(setdiff(c(1:5,7:9,12,15:20),c(2,5,7,12,15:20)),c(1:5,7:9,12,15:20),c(2,5,7,12,15:20),c(6,13,14), c(1:20), c(21:27)) # manual, indexes of genes
geneset_description_list = list('Vcor virulence genes (no toxins)','Vcor virulence genes','toxins','antibiotic resistance','all virulence genes','all zinc metalloproteases') # manual, for plot title

geneset_to_plot = 1 # manually change
gene_index_to_plot = gene_index_to_plot_list[[geneset_to_plot]]
geneset_description = geneset_description_list[[geneset_to_plot]]

## FOR LOOP: make a line plot for each gene of interest
for (z in 1:length(gene_list)){ # for each category
  
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
    
   #----- line + dot plot
   g = ggplot(data_plot, aes(x=time, y=log2fc)) + 
       geom_line(aes(color=gene)) + 
       geom_point(size=3,aes(color=gene, alpha=(padj<0.01))) + 
       geom_text(data=data_plot %>% filter(time == last(time)), aes(label = gene, 
                                                                  x = time - 15, 
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
        theme(legend.position='top')
       
    
  
    #----- SAVE figure (png first, for inspection)
    save_fn = paste0(save_folder,'/time_vs_log2fc_',geneset_description)
    # ggsave(paste0(save_fn,'.eps'), dpi=300, dev='eps', height = 3, width = 4, units = "in") 
    ggsave(paste0(save_fn,'.png'), dpi=300, dev='png', height = 6, width = 6, units = "in") 
    
} # end of for loop through each category of genes
