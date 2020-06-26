# Regurgitate DESeq2 results of single genes
# Input: JGI ID of gene of interest

# CREATED : 5/14/2020 by Cherry Gao



rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
load(file='dataset/jgi_downloaded/jgi_gene_data_clean.Rdata') # load cleaned & concatenated JGI gene ID data


# sequentially load and store log2FC 10m vs. 10c, and 60m vs. 60c DESeq2 results
# 10 min
load(file='deseq/deseq_kegg_concat/10muc_vs_10cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_10 = data_concat
rm(data_concat)
# 60 min
load(file='deseq/deseq_kegg_concat/60muc_vs_60cont_alpha_0.01_KEGGconcat.Rdata')
data_concat_60 = data_concat
rm(data_concat)

jgi_gene_of_interest = '647173052'

gene_data_10 = data_concat_10[which(data_concat_10$jgi_gene_id==jgi_gene_of_interest),]
gene_data_60 = data_concat_60[which(data_concat_60$jgi_gene_id==jgi_gene_of_interest),]

gene_data_10$padj
gene_data_60$padj

as.numeric(paste(gene_data_10$padj)) < 0.01
as.numeric(paste(gene_data_60$padj)) < 0.01

gene_data_10$foldchange
gene_data_60$foldchange

2^as.numeric(paste(gene_data_10$foldchange))
2^as.numeric(paste(gene_data_60$foldchange))