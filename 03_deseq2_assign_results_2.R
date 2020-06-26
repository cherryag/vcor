# Assign results of differential expression

# After running DESeq2 (deseq2_group_design.R): 
  # 1) run the 'results' function
  # 2) change FDR alpha from default alpha = 0.1
  # 3) MA-plot save
  # 4) concatenate with gene name and description, and save differentially expressed gene data as CSV

# CREATED : 1/9/2020 by Cherry Gao
# UPDATED : 3/13/2020 to v2 to change FDR alpha cutoff
#           future to do: automatically aggregate data for bar plot (next R script)
#           future to do: loop through each pairwise comparison
#           future to do: create top 50 up- and down-regulated genes for each pairwise comparison according to padj


# load DESeq2
# library('DESeq2') # run only first time

# set up
rm(list = ls()) # clear environment
alpha_change = 0.01 # FDR for 'results' function; default alpha = 0.1 (should change to FDR cut off you're planning to use)

# make / set saving directory
save_dir = paste0("deseq/deseq2_results_alpha_",alpha_change,"_logfc_cutoff_none/")
# dir.create(save_dir) # create folder (first time only)

#-----------------------------------


# the 2 conditions to compare [MANUALLY iterate through pairwise comparisons]
  # 0muc vs. 0cont
  # 10muc vs. 10cont
  # 60muc vs. 60cont
  # 10cont vs. 0cont
  # 60cont vs. 0cont
  # 60cont vs. 10cont
  # 10muc vs. 0muc
  # 60muc vs. 0muc
  # 60muc vs. 10muc
cond1 = "60muc" # "10muc" # "0muc" # "60muc"
cond2 = "10muc" # "60muc" # "0cont" # "60cont"

# set up environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')
# dir() # list all files in directory
load(file="deseq/dds_group_design.RData") # load dds and data


# assign results -- MANUAL ITERATION
# Indicate pairwise comparison terms in the terms after "group".
# The contrast argument of the function results takes a character vector of length three: 
# the name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator.
# The results function performs independent filtering by default using the mean of normalized counts as a filter statistic.
# Specify the 2 conditions to compare between as string, cond1 and cond2
res <- results(dds, contrast=c("group", cond1, cond2))

# print out results
head(res[order(res$padj),],10)
summary(res)
design(dds) # recall what the design formula used for analysis was
# summary(results(dds, contrast=c("group","10muc","10cont"),independentFiltering=FALSE)) # turn off independent filtering

# [optional but recommended] change FDR cutoff alpha (see independent filtering)
# Note: changing alpha doesn't change the number of genes in my dataset.
# The default significance level for independent filtering is 0.1; however, you should set this to the FDR cut off you are planning to use. 
# To change FDR cutoff alpha, independentFiltering=TRUE
# Changing alpha will change the results summary table (e.g. ratios).
res_adj <- results(dds, contrast=c("group", cond1, cond2), alpha=alpha_change) 
sum(res_adj$padj < alpha_change, na.rm=TRUE) # adjusted p values for the genes which do not pass the filter threshold (alpha) are set to NA.
sum(abs(res_adj$log2FoldChange) > 1 & res_adj$padj<alpha_change, na.rm = TRUE)

# save summary of results in a text file
sink(paste(save_dir,"summary_",cond1,"_vs_",cond2,"_alpha_",alpha_change,".txt", sep=""))
  summary(res_adj)
sink()


    # [optional] visualize independent filtering of results
    # The results function of the DESeq2 package performs independent filtering by default using the mean of normalized counts as a filter statistic. 
    # A threshold on the filter statistic is found which optimizes the number of padj lower than a specified significance level alpha.  
    # Genes with padj = na are ones DESeq2 has filtered out.
#    metadata(res05)$alpha # default alpha = 0.1 (should change to FDR cut off you're planning to use)
#    metadata(res05)$filterThreshold
#    plot(metadata(res05)$filterNumRej, 
#         type="b", ylab="number of rejections", 
#         xlab="quantiles of filter")
#    lines(metadata(res05)$lo.fit, col="red") 
#    abline(v=metadata(res05)$filterTheta)

# display metadata of results
mcols(res_adj,use.names=T) # function 'mcols' on the results object: meaning of the columns of res and what tests were used (res is a DataFrame object)

# MA-plot [save plot]
# Shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
# Points will be colored red if the adjusted p value is less than 0.1. 
# Points which fall out of the window (set by ylim) are plotted as open triangles pointing either up or down.
pdf(file=paste(save_dir,"MA_plot_", cond1, '_vs_', cond2, ".pdf",sep="")) # to save
plotMA(res_adj, ylim=c(-5,5), 
       main = paste("MA-plot", cond1, 'vs', cond2),
       xlab = "mean of normalized counts") 
# abline(h=c(-1,1),col="dodgerblue",lwd=2) # line at specific logFC
dev.off() # to save image

    # [optional] interactively label MA-plot
    # After calling plotMA, interactively detect row number of individual genes by clicking on individual data points on plot.
    # Recover the gene identifiers by saving the results indices.
#    idx <- identify(res05$baseMean, res05$log2FoldChange)
#    rownames(res05)[idx]

    # [optional] plot single gene cross groups
    # function plotCounts normalizes counts by sequencing depth and adds a pseudocount of 1/2 to allow for log scale plotting.
    # The counts are grouped by the variables in intgroup, where more than one variable can be specified.
    # Here we specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by rowname or by numeric index.
#    plotCounts(dds, gene=which.min(res05$padj), intgroup="condition")

# merge results with gene names + export as csv
# rm(degenelist,genenames,rownames) # overwrite on this variable
genenames <- counts_raw_original[,1:3] # geneID, geneLength, and geneProduct columns
  rownames <- res_adj@rownames # geneIDs
  baseMean <- res_adj$baseMean
  foldchange <- res_adj$log2FoldChange
  lfcSE <- res_adj$lfcSE
  stat <- res_adj$stat
  pvalue <- res_adj$pvalue
  padj <- res_adj$padj
padj_table <- cbind(rownames,baseMean,foldchange,lfcSE,stat,pvalue,padj) # correct formatting for merging
merged_table <- merge(padj_table,genenames,by.x="rownames",by.y="geneID")
degenelist <- merged_table[order(as.numeric(as.character(merged_table$padj))),] # order by padj (convert from 'factor' to 'numeric')


# [optional] append gene product name ot res_adj
# load cleaned Vcor genome data downloaded from JGI + extract first element of each list in data_jgi (gene ID)
load(file = 'dataset/jgi_downloaded/jgi_gene_data_clean.Rdata')
jgi_gene_list = unlist(lapply(data_jgi,'[[',1)); 

# loop for append gene product with res_adj
gene_descriptor_list_temp = matrix(ncol = 1, nrow = length(res_adj$padj)) # initialize temp container
for (j in 1:length(res_adj$padj)) { 
  line = res_adj[j,] # process each line
  gene_oid = rownames(line) # extract JGI gene ID
  # extract JGI gene ID descriptor 
  if (gene_oid %in% jgi_gene_list) {
    i = match(gene_oid,jgi_gene_list) # extract index
    if (gene_oid != data_jgi[[i]][1]){'expression data gene ID is not equal to JGI gene ID!'} # sanity check
    gene_descriptor = data_jgi[[i]][5] # 5th column = gene descriptor
    #   gene_descriptor_list_temp[i] <- gene_descriptor # append to temp container
    gene_descriptor_list_temp[j] <- paste(gene_descriptor,gene_oid,sep='_')
  }
} # end of cycling through line numbers
# append descriptor as a new column
res_adj$product <- gene_descriptor_list_temp


# write csv
write.csv(as.data.frame(degenelist),
          file = paste(save_dir,cond1,"_vs_",cond2,"_alpha_",alpha_change,".csv", sep=""))

# save .Rdata [OPTIONAL]
save(degenelist,res_adj, file = paste(save_dir,cond1,"_vs_",cond2,"_FDRalpha_",alpha_change,".Rdata", sep=""))

