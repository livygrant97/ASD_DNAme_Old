################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

packages <- c("missMethyl", "minfi", "cluster","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","dplyr","bigmelon","clusterProfiler","DMRcate","ggplot2","factoextra","tidyverse")
lapply(packages, library, character.only = TRUE)



###############################################################
#### loading and cleaning the data ####
###############################################################

#read in results
if(file.exists("/path/to/file/dasenAutosomeCleaned.Rata")){
  load("/path/to/file/dasenAutosomeCleaned.Rata")
} else{
diffMF<-read.csv("/path/to/file/diff_cpg_between_F_M.xls")
#performing FDR test
diffMF$FDR<-p.adjust(diffMF$P.Value,method="fdr")
#subset those probes which FDR <0.05
Passed_FDR<- subset(diffMF,FDR < 0.05)
#load annotations
data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)
#get probes on x and y chromosome
sexProbes<-as.character(annotationTable$Name[annotationTable$chr %in% c("chrX","chrY")])
#remove probes not on autosomes
sigCpG_Autosomes<- Passed_FDR[!Passed_FDR$Row.names %in% sexProbes,]
#read in annotations for SNP related probes remove them from analysis
SNPAnnotations<-read.delim("/home/og16379/diff_cpg_fm/data/EPIC.hg19.manifest.tsv")
SNPAssociatedCpG<- SNPAnnotations %>% dplyr::filter(MASK_general==TRUE)
dasenAutosomeCleaned<- sigCpG_Autosomes[!sigCpG_Autosomes$Row.names %in%  SNPAssociatedCpG$probeID,]
# downloaded the list of cross hybridizing probes from here https://pubmed.ncbi.nlm.nih.gov/23314698/
chProbes<-read.table("path/to/file.txt")
dasenAutosomeCleaned<-dasenAutosomeCleaned[!(dasenAutosomeCleaned$Row.names
  %in% chProbes$V1),]
paste(formatC(nrow(SNPAssociatedCpG),big.mark = ','),'Number of SNP associated probes or cross-hybridizing probes removed')
save(dasenAutosomeCleaned,file="/path/to/file/dasenAutosomeCleaned.Rata")
}



#################################################################################
#### Analysing the data ####
#################################################################################


#colours we use for plotting throughout the manuscript
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
maleColour<-cbbPalette[7]
femaleColour<-cbbPalette[6]

###############################################################
####  Generate Table 1 ####
###############################################################

# character vector of all CpG sites tested.
all.cpg <- as.character(rownames(annotationTable))
# character vector of significant CpG sites to test for GO term enrichment
sig.cpg<-as.character(dasenAutosomeCleaned$Row.names))
# perform gene ontology testing
sexAssociatedSitesPathways <- gometh(sig.cpg,all.cpg,collection="GO",array.type="EPIC")
#get signficiant pathways
table(sexAssociatedSitesPathways$FDR<0.05)
#table of top GO results
topGO(sexAssociatedSitesPathways)


###############################################################
####  Filtering list to get strongest differences ####
###############################################################

# then we wanted to filter this list down to select the CpGs which showed the strongest differences
# between males and females
# first, detirmine what sex methylation is higher in for each probe using the average beta values
# add a column called $Higher_Met_In which indicates which sex methylation is higher in
dasenAutosomeCleaned$Higher_Met_In<-"Male"
dasenAutosomeCleaned$Higher_Met_In [dasenAutosomeCleaned$F_AVG >dasenAutosomeCleaned$M_AVG] <- "Female"
#split for sex so you have all probes regardless of deltaBeta. incase want to change threshold
splitCpG_sig<-split(dasenAutosomeCleaned,dasenAutosomeCleaned$Higher_Met_In)
femalesSigCpG<-splitCpG_sig$Female
malesSigCpG<-splitCpG_sig$Male
# we used threshold of delta beta diff greater than 0.05
# make object for each for further analysis
sigProbesFemales<- subset(femalesSigCpG, deltaBeta > 0.05  | deltaBeta < -0.05)
sigProbesMales<- subset(malesSigCpG, deltaBeta > 0.05 | deltaBeta <-0.05)
allSigProbes<-rbind(sigProbesMales,sigProbesFemales)
