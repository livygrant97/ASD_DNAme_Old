################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#below is the code to calculate DEG between males and females
#count matrix already pre-processed

################################################################################
# libraries
################################################################################

if(!require("DESeq2", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("DESeq2")
}
library(DESeq2)

# read in data
countMatrix<-read.delim("/home/og16379/diff_cpg_fm/data/GSE120312_Counts_Matrix.txt")
# add in identifiers for sex, confirmed with PCA
countData <- data.frame(row.names=(countMatrix$Ensembl.Gene.ID),countMatrix[,11:30])
sex<-c("Female","Male","Male","Male","Female","Female","Male","Male","Female","Male","Female","Female","Female","Female",
"Male","Male","Female","Male","Female","Male")
# construct meta data
metaData <- data.frame("id"=colnames(countData),"dex"=sex)

#construct deseqdata obj
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design=~dex)
#calculate DEG
dds <- DESeq(dds)
res <- results(dds)
head(results(dds))
summary(res)
res <- res[order(res$padj),]
head(res)
