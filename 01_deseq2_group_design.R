# Setup and run DESeq2
# "~ group" design formula instead of using interaction terms.
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# 1) load data
# 2) set up design formula 
# 3) generate variable 'dds'
# 4) save workspace environment

# CREATED : 1/3/2020 by Cherry Gao


rm(list = ls()) # clear environment
getwd() # current directory

# set directory
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal/')
dir() # list all files in directory

# installation (first time only)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

# load DESeq2
library("DESeq2")

# load data (raw counts data from JGI - should be unnormalized)
  # colData was made by hand (maps sample name to muc/cont and time point)
  # countData has sample names (columns) and geneID (rows) and the raw count data
countData_raw = read.table('countData.txt', sep='\t', header=T, row.names=1)
colData_raw = read.table('colData.txt', sep='\t', header=T, row.names=1)
counts_raw_original = read.delim("raw_counts_JGI.txt") # table including gene names; something goes wrong here if use read.table (import dataset using UI instead)

# check that sample names are in the same order (i.e. colname of countData and rowname of colData)
all(rownames(colData_raw) %in% colnames(countData_raw))
all(rownames(colData_raw) == colnames(countData_raw))

# choose a reference level - control group should be set as the first "factor level"
  # Tell DESeq2 functions which level represents the control group (by default, R will choose reference levels for factors based on alphabetical order).
  # By setting the control as first factor level, default log2 fold changes are calculated as muc over control.
  # Also make time a factor s.t. model doesn't assume numeric values have increasing fold change for higher values 
  # Alternatively, use the "contrast" argument in results.
colData_raw$time <- factor(colData_raw$time)
colData_raw$condition <- factor(colData_raw$condition, levels=c('cont','muc'))

# construct DESeqDataSet -- group design formula instead of using interaction terms
  # Pay special attention to design formula (starting with ~).
  # Design formula expresses the variables used for modeling.
  # "design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model."
  # section "Interactions"
    # 1) "combine the factors of interest into a single factor with all combinations of the original factors"
    # 2) "change the design to include just this factor, e.g. ~ group"
  # "Using this design ("~ group") is similar to adding an interaction term, in that it models multiple condition effects which can be easily extracted with results."
  # In factor colData, levels are cont0, cont10, cont60, muc0, muc10, muc60
dds <- DESeqDataSetFromMatrix(countData = countData_raw,
                              colData = colData_raw,
                              design = ~ time + condition + condition:time)
#                              design = ~ condition + time + condition:time) # gives same result as first
#                              design = ~ time + condition + time:condition) # gives same result as first

dds$group <- factor(paste0(dds$time, dds$condition)) 
design(dds) <- ~ group

# perform differential expression analysis
  # dds = "DESeqDataSet" object that contains associated design formula
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# save
save.image(file='dds_group_design.RData')

# NEXT: deseq2_assign_results.R
