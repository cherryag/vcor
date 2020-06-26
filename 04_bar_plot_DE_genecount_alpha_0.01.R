# bar plot: counts of significantly differentially expressed genes 
# data is manually entered

# 3/14/2020 : manual entry for FDR alpha = 0.01


# [MANUAL ENTRY] number of genes LFC > 0 (up) and LFC < 0 (down), read from summary.txt files
up0 = 0 # 0M vs 0C
down0 = 0
up10 =  1379 # 10M vs 10C
down10 = 1326
up60 =  1159 # 60M vs 60C
down60 = 1076

# Pairwise controls
up_c10_c0 = 797 # 10C vs. 0C
down_c10_c0 = 847
up_c60_c0  = 224 # 60C vs. 0C
down_c60_c0 = 285
up_c60_c10 = 200 # 60C vs. 10C
down_c60_c10 = 299 

# Pairwise mucus
up_m10_m0 = 1599 # 10M vs. 0M
down_m10_m0 = 1486
up_60m_m0 = 1324 # 60M vs. 0M
down_60m_m0 = 1374 
up_m60_m10 = 1307 # 60M vs. 10M
down_m60_m10 = 1334


# aggregate gene counts into labeled matrix
up <- c(up0,up10,up60,   up_c10_c0,up_c60_c0,up_c60_c10,    up_m10_m0,up_60m_m0,up_m60_m10)
down <- c(down0,down10,down60,   down_c10_c0,down_c60_c0,down_c60_c10,    down_m10_m0,down_60m_m0,down_m60_m10)
counts <- matrix(c(down,up),nrow = length(up),ncol = 2)
counts <- t(counts) # transpose
rownames(counts) <- c("downregulated","upregulated")

# plot
pdf(file=paste(save_dir,"bar_DE_genecounts_all.pdf",sep=""), width=6, height=5) # to save -- need to be polished
        barplot(counts, 
                main=paste("Number of significantly differentially expressed genes","\n(log2 fc threshold = none; padj < )",alpha_change),
                names.arg = c("im / ic","10m / 10c","60m / 60c",   "10c / ic","60c / ic","60c / 10c",     "10m / im","60m / im","60m / 10m"),
                ylab = "Number of genes",
                legend.text = c(rownames(counts)),
                args.legend = list(x = "top" , bty = "n"),
                cex.names=0.75, las=2)
dev.off() # to save image