##Heatmap plotting
install.packages("reshape")
install.packages("ggplot2")

setwd("~/Documents/Spaceplankton/Transcriptomics/DE_genes/")

library(reshape)
library(ggplot2)
library(scales)

DataTable <- read.csv("Transcripts_of_interest_for_paper_for_plotting.txt", sep="\t", header=TRUE)
PlotTable <- DataTable[which(DataTable$fdr<=0.5),c(1,3,7:21)]
PlotTable <- within(PlotTable, Enzyme <- paste(ID,enzyme, sep="_"))
PlotTable$Enzyme <- factor(PlotTable$Enzyme, levels=PlotTable$Enzyme)
PlotTable <- PlotTable[,-c(1:2)]
#normalize values to mean of set1 ground samples (new values are fold-changes), 
#and take the log2 to plot log2(fold-changes)
for (i in 1:length(row.names(PlotTable))) {
  PlotTable[i,3:14] <- log(PlotTable[i,3:14]/mean(as.matrix(PlotTable[i,9:11])),2)
  #PlotTable[i,3:14] <- PlotTable[i,3:14]/max(PlotTable[i,3:14])
  #divAVGtable$max[i] <- max(as.matrix(divAVGtable[i,3:14]))
}
tidy_table <- melt(PlotTable, id=c("Enzyme","function.","fdr","fdr.set2"))

#create vector with row numbers of all entries with fdr<0.1 AND only those functions to plot
significantlyDEstress1    <- which(PlotTable[grep("A|B|C|D|F|G|H|I", PlotTable$function.),]$fdr<=0.1)
significantlyDEflagellum1 <- which(PlotTable[grep("J|K", PlotTable$function.),]$fdr<=0.1)
significantlyDEsec.met1   <- which(PlotTable[grep("L_", PlotTable$function.),]$fdr<=0.1)
significantlyDEstress2    <- which(PlotTable[grep("A|B|C|D|F|G|H|I", PlotTable$function.),]$fdr.set2<=0.1)
significantlyDEflagellum2 <- which(PlotTable[grep("J|K", PlotTable$function.),]$fdr.set2<=0.1)
significantlyDEsec.met2   <- which(PlotTable[grep("L_", PlotTable$function.),]$fdr.set2<=0.1)
#create vectors with row numbers of only one set (1 or 2) AND only those functions to plot
stress1    <- intersect(grep("set1", tidy_table$variable),grep("A|B|C|D|F|G|H|I", tidy_table$function.))
flagellum1 <- intersect(grep("set1", tidy_table$variable),grep("J|K", tidy_table$function.))
sec.met1   <- intersect(grep("set1", tidy_table$variable),grep("L_", tidy_table$function.))
stress2    <- intersect(grep("set2", tidy_table$variable),grep("A|B|C|D|E|F|G|H|I", tidy_table$function.))
flagellum2 <- intersect(grep("set2", tidy_table$variable),grep("J|K", tidy_table$function.))
sec.met2   <- intersect(grep("set2", tidy_table$variable),grep("L_", tidy_table$function.))

#calculate the average of the Set1 ground samples to set the color scale midpoint to this average


#quick heatmap (to plot the different heatmaps change dataframe tidy_table[...,], the min/max values,
# as well as the annotation data x=...)
breaks_stress=c(-4,-3,-2,-1,0,1,2,3,4)
n <- length(tidy_table$value[sec.met1])
secondhighest_value<-sort(tidy_table$value[sec.met1])[n-1]
secondlowest_value<-sort(tidy_table$value[sec.met1],decreasing=TRUE)[n-1]

#definition of the color scale scaled on the min or max log2FC, depending on which is closer to zero. 
#keep color scale symmetrical, saturating at the abs value closest to 0+0.2
myplot <- ggplot(tidy_table[sec.met1,], aes(Enzyme, variable, fill = value)) + geom_tile() +
  scale_fill_gradientn(colors=c("blue","white","red"),
            values=rescale(c(min(tidy_table$value[sec.met1]),0,max(tidy_table$value[sec.met1]))),
            limits=c(min(tidy_table$value[sec.met1]),max(tidy_table$value[sec.met1])),breaks=breaks_stress) +  
  #midpoint jetzt auf log(1)=0, weil Daten logarithmiert sind und normalisiert auf mean(ground_set1)
  scale_x_discrete(position = "top") + 
  #geom_vline(xintercept = c(11.5,28.5,47.5,54.5),color = "black", size=0.3) +
  #geom_vline(xintercept = 11.5,color = "black", size=0.5) +
  geom_hline(yintercept = 3.5,color = "black", size=0.3,linetype = "longdash") +
  annotate("text", label=" *",x=significantlyDEsec.met1-0.3,y=0,size=5)
myplot + coord_flip() +
  theme(axis.text.x.bottom =element_text(size=10,angle=270, hjust = 0, vjust=0.3,colour="grey30")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+ labs(fill="log2(fold-change)")


pdf(file="..\\Figures_Results\\secmet1_map.pdf", 7, 3)
myplot + coord_flip() +
  theme(axis.text.x.bottom =element_text(size=10,angle=270, hjust = 0, vjust=0.3,colour="grey30")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
dev.off()

#note meeting Jeanette: can also use a saturating color scale instead of taking the log (like full color at 1, everything above is full color)

