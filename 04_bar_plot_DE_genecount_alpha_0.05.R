# bar plot: counts of significantly differentially expressed genes 

# [MANUAL ENTRY] number of genes LFC > 0 (up) and LFC < 0 (down), read from summary.txt files
up0 = 0 # TIM vs TIC
down0 = 0
up10 = 1570 # T10M vs T10C
down10 = 1545
up60 = 1386 # T60M vs T60C
down60 = 1319

up_c0_c10 = 1100 # TIC vs. T10C
down_c0_c10 = 1063
up_m10_m60 = 1499 # T10M vs. T60M
down_m10_m60 = 1530

# aggregate gene counts into labeled matrix
up <- c(up0,up10,up60,up_c0_c10,up_m10_m60)
down <- c(down0,down10,down60,down_c0_c10,down_m10_m60)
counts <- matrix(c(down,up),nrow = length(up),ncol = 2)
counts <- t(counts) # transpose
rownames(counts) <- c("downregulated","upregulated")

# plot
pdf(file=paste("bar_DE_genecounts_all.pdf",sep=""), width=6, height=5) # to save -- need to be polished
        barplot(counts, 
                main=paste("Number of significantly differentially expressed genes","\n(log2 fc threshold = none; padj < 0.05)"),
                names.arg = c("TIM vs TIC","T10M vs T10C","T60M vs T60C","TIC vs T10C","T10M vs T60M"),
                ylab = "Number of genes",
                legend.text = c(rownames(counts)),
                args.legend = list(x = "topleft" , bty = "n"),
                cex.names=0.75, las=2)
dev.off() # to save image