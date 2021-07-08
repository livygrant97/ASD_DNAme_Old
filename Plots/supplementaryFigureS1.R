################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

# below is all the code to produce the supplementary figure S1
# previous code must have been run in order to produce these figures

################################################################################
# Generate figure S1A
################################################################################

 pvals<- sigCpG_Autosomes_Cleaned$p
 
 dmpQQ <- function(pvals){

            fold.range <- c(0,6)
             prob.range <- c(0,300)
             plot(
                     x = 0, y = 0, xlim = fold.range, ylim = prob.range,
                     xlab = "", ylab = "",
                     axes = TRUE, type = "n", cex.axis = 1.0, cex.lab = 1.0)

             coords <- par("usr")
             gx <- grconvertX(coords[1:2], "user", "in")
             gy <- grconvertY(coords[3:4], "user", "in")
             width <- max(gx) - min(gx)
             height <- max(gy) - min(gy)

             tmp <- tempfile()
             png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")

             plot.new()
             par(mar = c(0.0, 0.0, 0.0, 0.0))
             qq(pvals)
             dev.off()
             panel <- readPNG(tmp)
             rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
             mtext(text = "Expected -log10 (p)", side = 1, line = 4, cex = 1)
             mtext(text = "Observed-log10(p)", side = 2, line = 2.5, cex = 1)

    }

   pdf("/home/og16379/diff_cpg_fm/pdf/rasterisedQQ.pdf",width=8)
     dmpQQ(pvals)
 dev.off()


 ################################################################################
 # Generate figure S1B and perform wilcoxon test
 ################################################################################

 #first estimate cell type propotions
 ecc<-estimateCellCounts.gds(gfile,referencePlatform="IlluminaHumanMethylationELoading required package: FlowSorted.Blood.450k"))

 #remove outliers
 eccDf<-as.data.frame(ecc)
 eccDf<-eccDf[info$barcode,]
 eccDf$Sex<-info$nsex
 eccDf$Sex<- gsub('F', 'Females', eccDf$Sex)
 eccDf$Sex<- gsub('M ', 'Males', eccDf$Sex)

 eccMelt<-melt(eccDf)

 #plot
 ggplot(eccMelt,aes(x=variable,y=value,fill=Sex))+
 geom_boxplot()+
 labs(x="Blood cell types",y="Proportion")+theme(legend.title=element_blank())+
   theme(axis.line = element_line(colour = "black"),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.background=element_rect(fill="white",colour="white"),
     panel.border = element_rect(fill=NA,color="black"))+
     scale_fill_manual(values=c(femaleColour,maleColour))

 #wilcoxon p values
 f<-split(eccMelt,eccMelt$Sex)
 femaleData<-f$Females
 maleData<-f$Males
 maleData<-split(maleData,maleData$variable)
 femaleData<-split(femaleData,femaleData$variable)

 maleMatrix<-matrix(0,486,6)
 colnames(maleMatrix)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
 maleMatrix[,1]<-maleData$CD8T$value
 maleMatrix[,2]<-maleData$CD4T$value
 maleMatrix[,3]<-maleData$NK$value
 maleMatrix[,4]<-maleData$Bcell$value
 maleMatrix[,5]<-maleData$Mono$value
 maleMatrix[,6]<-maleData$Gran$value

 femaleMatrix<-matrix(0,685,6)
 colnames(femaleMatrix)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
 femaleMatrix[,1]<-femaleData$CD8T$value
 femaleMatrix[,2]<-femaleData$CD4T$value
 femaleMatrix[,3]<-femaleData$NK$value
 femaleMatrix[,4]<-femaleData$Bcell$value
 femaleMatrix[,5]<-femaleData$Mono$value
 femaleMatrix[,6]<-femaleData$Gran$value

 #repeat in loop
 wilcox.test(femaleMatrix[,1],maleMatrix[,1])
