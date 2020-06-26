# Heatmap of log2FC gene expression


#---------- HISTORY ----------
# CREATED : 5/13/2020 (modified from Uli's heatmap script)
# UPDATED : 7/16/2020 (added 3rd option for data input, directly from KEGG pathway set Excel sheet)
#-----------------------------

#----------SUMMARY----------
# 1) organize data for plotting
# 2) plot
#---------------------------

# install packages
# install.packages("reshape")
# install.packages("ggplot2")

library(reshape)
library(ggplot2)
library(scales)
library(readxl)

#---------- SET UP
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
save_dir='/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/virulence_genes/heatmap_log2fc' # where to save timelapse plots
  # dir.create(save_dir)

#------ INPUT DATA OPTION 1: .Rdata
# load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data
# load(file='dataset/virulence_genes_curated/jgi_gene_id_vps.Rdata') # VPS biofilm genes

# gene_list = jgi_id_and_name_manual
# geneset_description = 'secretion systems' # manual
#----------------------------------


#------- INPUT DATA OPTION 2: Excel
# geneset_description = 'coralliilyticus pathogenicity islands' # [5/15/2020]
# geneset_description = 'gene_list_chemotaxis' # [5/15/2020]
# geneset_description = 'gene_list_MCP' # [5/15/2020]
# geneset_description = 'gene_list_flagella' # [5/15/2020 + 5/28]
# geneset_description = 'vps' # [5/17/2020]
# geneset_description = 'gene_list_secretion_system' 
# geneset_description = 'virulence_genes_ab_resistance'# [5/28/2020]
# geneset_description = 'gene_list_nqr' # [5/28/2020]
# geneset_description = 'gene_list_flagella_che' # [5/28/2020]

# v = read_excel('dataset/virulence_genes_curated/secretion_system_genes.xlsx') # [ADDED 5/2/2020] circumvent having to process in virulence_genes_manually_curated.R
  # geneset_description = 'secretion systems' # manual
# v = read_excel(paste0('dataset/virulence_genes_curated/',geneset_description,'.xlsx')) 

# assign gene_list (JGI + gene name)
# gene_list = cbind(v$`JGI #`,v$`gene name`)
#----------------------------------

#------- INPUT DATA OPTION 3: KEGG Pathway list of JGI ID #
load('dataset/kegg_downloaded/kegg_pathways_jgi.Rdata')
kegg_path_num = list('03010','00970') # MANUAL
geneset_description ='ribosome_trna_synthesis' # MANUAL

# assign gene_list (JGI ID only)
gene_list = list()
# for loop to construct list of JGI ID numbers from different pathways
for (i in 1:length(kegg_path_num)){ 
  path_to_append = kegg_list_jgi[[kegg_path_num[[i]]]]
  gene_list = append(gene_list, path_to_append[2:length(path_to_append)]) # eliminate the first row 
}
gene_list = t(t(as.numeric(unlist(gene_list)))) # convert to numeric, then transpose twice (to be able to index with [,1])
#---------------------------------

# clean gene_list
index_to_eliminate = matrix() # initialize the list
index_to_eliminate = which(is.na(gene_list[,1]))
index_to_eliminate = c(index_to_eliminate,  which(gene_list[,1] == '-'))
if (length(index_to_eliminate)>0){ gene_list = gene_list[-index_to_eliminate,] } # eliminated genes with no JGI gene ID assignment

# sequentially load and store log2FC 10m vs. 10c, and 60m vs. 60c DESeq2 results
# 10 min
load(file='deseq/deseq_kegg_concat/10muc_vs_10cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_10 = data_concat
rm(data_concat)
# 60 min
load(file='deseq/deseq_kegg_concat/60muc_vs_60cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_60 = data_concat
rm(data_concat)

# extract indexes of genes in data_concat 
gene_idx_10 = match(as.numeric(unlist(gene_list[,1])), data_concat_10$seq_data_gene_id)
gene_idx_60 = match(as.numeric(unlist(gene_list[,1])), data_concat_60$seq_data_gene_id)

# append gene product name from data_concat if INPUT DATA OPTION 3
if (ncol(gene_list) == 1){
  gene_list = as.data.frame(gene_list)
  gene_list$'product' = data_concat_10$Product_name[gene_idx_10]}

# order of plotted bars (gene_order = index of gene_idx_10)
gene_order = c(1:nrow(gene_list)) # arbitrarily linear
# gene_order = order(as.numeric(paste(data_concat_10$foldchange[gene_idx_10]))) # order by log2fold change value (for MCP figure)

gene_idx_10_ordered = gene_idx_10[gene_order]
gene_idx_60_ordered = gene_idx_60[gene_order]

# gene names - manual (in excel) or from database
# gene_name_list = paste0(data_concat_10$Locus_tag[gene_idx_10_ordered],'_',data_concat_10$Product_name[gene_idx_10_ordered]) # from JGI database
gene_name_list = paste0(data_concat_10$Locus_tag[gene_idx_10_ordered],'_',gene_list[,2]) # from manual excel sheet

# concatenate gene data to plot
fc_10 = matrix(); padj_10 = matrix(); fc_60 = matrix(); padj_60 = matrix(); data_plot_temp = matrix() # initialize

fc_10 = as.numeric(paste(data_concat_10$foldchange[gene_idx_10_ordered]))
padj_10 = as.numeric(paste(data_concat_10$padj[gene_idx_10_ordered]))
  data_10 = cbind(10, fc_10, padj_10, as.character(gene_name_list), c(1:length(gene_name_list)))

fc_60 = as.numeric(paste(data_concat_60$foldchange[gene_idx_60_ordered]))
padj_60 = as.numeric(paste(data_concat_60$padj[gene_idx_60_ordered]))
  data_60 = cbind(60, fc_60, padj_60, as.character(gene_name_list), c(1:length(gene_name_list)))

# format gene data for plotting
data_plot_temp = rbind(data_10,data_60)
  
data_plot <- data.frame(
  time = as.numeric(data_plot_temp[,1]),
  log2fc = as.numeric(data_plot_temp[,2]),
  padj = as.numeric(data_plot_temp[,3]),
  gene = data_plot_temp[,4],
  order = as.numeric(data_plot_temp[,5])
)

# extract indexes of significant and non-significant DE genes
alpha = 0.05 # cutoff
s_index = which(data_plot$padj < alpha) # significant 
ns_index = which(data_plot$padj >= alpha) # not significant


# heatmap
g <- ggplot(data_plot, aes(x=order, y=time, fill=log2fc)) + geom_tile() + #scale_x_reverse() +
    scale_x_reverse(breaks = 1:length(gene_name_list), labels = data_plot$gene[1:length(gene_name_list)]) + # label genes, order linearly
    # scale_x_reverse(breaks = 1:length(gene_name_list), labels = data_plot$gene[gene_order]) + # label genes, order by gene_order
# g <- ggplot(data_plot, aes(gene, time, fill = log2fc)) + geom_tile() +
     scale_fill_gradientn(colors=c("purple","white","orange"), # c("#bdbdbd", "#f7f7f7", "#f8992e") / c("#f1a340", "#f7f7f7", "#998ec3") # purple/white/orange from colorbrewer

# COLORBAR, BY MIN AND MAX VALUES:
#     values=rescale(c(min(data_plot$log2fc),0,max(data_plot$log2fc))), # scale with 0 in the middle
#     limits=c(min(data_plot$log2fc),max(data_plot$log2fc)))  

# COLORBAR, SYMMETRICAL BY ABSOLUTE VALUE OF MIN OR MAX
  values=rescale(c(-max(data_plot$log2fc),0,abs(max(data_plot$log2fc)))), # scale with 0 in the middle
  limits=c(-max(data_plot$log2fc),abs(max(data_plot$log2fc))))  
      
# COLORBAR, MANUAL RANGE
 #   values=rescale(c(-7.2,0,5.5)), limits=c(-7.2,5.5))  

# mark significant genes with '*' and non-significant genes with 'x'
g <- g + annotate("text", label=" *", x = data_plot$order[s_index], y=data_plot$time[s_index], size = 5) # + 
       #  annotate("text", label=" x", x = data_plot$order[ns_index], y=data_plot$time[ns_index], size = 4)

#---------- OPTIONAL cosmetic things ----------
  
  # change the side of the gene label
  # g <- g + scale_x_discrete(position = "top") # change the side of the gene label

  # horizontal line to divide genes
  # g <- g + geom_vline(xintercept = c(41)+0.5, color = "black", size=0.3) # c(1,12,21,29,37,49)+0.5 secretion system

  # dashed line between the columns
  g <- g + geom_hline(yintercept = 35, color = "black", size=0.3, linetype = "longdash")
   
#---------------------------------------------

g <- g + coord_flip() 

# get rid of grey background + colorbar label + title
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),                      # get rid of grey background
               panel.background = element_blank(), axis.line = element_blank()) + labs(fill="log2fc") +     # colorbar label
         labs(title = paste('Heatmap: log2 fold change of ',geneset_description, 'genes'),
              subtitle = paste('x = not significant (padj ≥', alpha, '); * = significant \n color scaled on all (n.s. and s.)'))
     


#----- SAVE figure
save_fn = paste0(save_dir,'/heatmap_log2fc_',geneset_description)
ggsave(paste0(save_fn,'_same_cbar_large.pdf'), dpi=300, dev='pdf', height = 30, width = 30, units = "in", useDingbats=FALSE) # useDingbats for ggplot -> Adobe Illustrator transition (markers are not fonts)
ggsave(paste0(save_fn,'_same_cbar_large.png'), dpi=300, dev='png', height = 30, width = 30, units = "in") 

