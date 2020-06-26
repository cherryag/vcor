# Data quality assessment by sample clustering and visualization
# After running DESeq2, make heat maps of sample-to-sample distances and PCA. 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  # "In order to test for differential expression, we operate on raw counts and use discrete distributions 
  # as described in the previous section on differential expression. However for other downstream analyses 
  # - e.g. for visualization or clustering – it might be useful to work with transformed versions of the 
  # count data."

##### installation (first time only) ######
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("pheatmap")
###########################################

# set up environment
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/')
load(file="deseq/dds_group_design.RData") # load dds and data

# make / set saving directory
save_dir = paste0("deseq/overall_data_vis/")
# dir.create(save_dir) # create folder (first time only)

# 1) log2(n + 1)
ntd <- normTransform(dds)

# Two other normalization methods: 
  # "Both transformations produce transformed data on the log2 scale which has been normalized with 
  # respect to library size or other normalization factors."

  # "The point of these two transformations, the VST and the rlog, is to remove the dependence of the 
  # variance on the mean, particularly the high variance of the logarithm of count data when the mean is low.
  # Both VST and rlog use the experiment-wide trend of variance over mean, in order to transform the data to remove 
  # the experiment-wide trend. Note that we do not require or desire that all the genes have exactly the same 
  # variance after transformation. Indeed, in a figure below, you will see that after the transformations the genes 
  # with the same mean do not have exactly the same standard deviations, but that the experiment-wide trend has flattened. 
  # It is those genes with row variance above the trend which will allow us to cluster samples into interesting groups."

# 2) variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010)
# Note: vst with blind=FALSE was recommended on a forum by Michael Love
vsd <- vst(dds, blind=FALSE)

# 3) regularized logarithm (rlog), which incorporates a prior on the sample differences (Love, Huber, and Anders 2014)
rld <- rlog(dds, blind=FALSE)

# indicate normalization method for figure generation
data <- vsd
norm_method = "variance stabilizing transformation" # "log2(n+1)" # for titles
norm_method_savename = "vst"

#---------- heatmap of the count matrix (top 20 genes with most counts)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","time")])

# gene heatmap figure
pdf(file = paste0(save_dir,"heatmap_counts_samples_",norm_method_savename,".pdf"), width=8, height=5) # to save
  pheatmap(assay(data)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=df,
           main=paste("heatmap of count matrix of top 20 genes with the most counts","\nsummed across all samples (normalization method: ",norm_method, ")"))
dev.off()

# save gene heatmap data
write.csv(as.data.frame(assay(data)[select,]),
          file = paste0(save_dir,"heatmap_counts_samples_",norm_method_savename,".csv"))

  

#---------- heatmap of sample-to-sample distances
library("RColorBrewer")
sampleDists <- dist(t(assay(data))) # dist function to the transpose of the transformed count matrix to get sample-to-sample distances; assay function is used to extract the matrix of normalized values
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(data@colData)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# save sample distance heatmap data
write.csv(as.data.frame(sampleDistMatrix),
          file = paste0(save_dir,"sample_distance_heatmap_",norm_method_savename,".csv"))

# sample distance heatmap figure
pdf(file = paste0(save_dir,"sample_distance_heatmap_",norm_method_savename,".pdf"), width=5, height=5) # to save -- need to be polished
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           main=paste("sample-to-sample distances", "\nnormalization method: ",norm_method))
dev.off() # to save image



#---------- PCA of the samples - customize using ggplot function (first return data)
library("ggplot2")
pcaData <- plotPCA(data, intgroup=c("condition", "time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 1)

# save PCA data
write.csv(as.data.frame(pcaData),
          file = paste0(save_dir,"pca_data_",norm_method_savename,".csv"))

# PCA plot
pdf(file=paste0(save_dir,"pca_",norm_method_savename,".pdf"), width=5, height=5) # to save
  z <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  nudge <- position_nudge(x = 5, y = 0) # 'nudge' the points up a bit
  z + 
    geom_text(aes(label = name), position = nudge, size = 1.2) + # geom_label = label dots with box; geom_text = label without box
    theme_bw() + # get rid of background 
    labs(title=paste("PCA: ",norm_method)) # title
dev.off()
