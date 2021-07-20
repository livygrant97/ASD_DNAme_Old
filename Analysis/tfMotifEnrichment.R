################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#below is the code to perform the tf motif enrichment

################################################################################
# libraries
################################################################################
if(!require("MotifDb", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("MotifDb")
}
library(MotifDb)


if(!require("PWMEnrich", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich")
}
library(PWMEnrich)



if(!require("PWMEnrich.Hsapiens.background", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich.Hsapiens.background")
}
library(PWMEnrich.Hsapiens.background)


if(!require("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)


################################################################################
# get data
################################################################################
# load in data
data(PWMEnrich.Hsapiens.background)
data(PWMLogn.hg19.MotifDb.Hsap)
# register cores
registerCoresPWMEnrich(30)

# load in bed files for the saDMPs as GRanges
femaleGRanges<-import("/home/og16379/diff_cpg_fm/data/resFemale.bed")
maleGRanges<-import("/home/og16379/diff_cpg_fm/data/resMale.bed")

# extend the GRanges by 100bp
femaleGRanges<-resize(femaleGRanges,fix="start",width=width(femaleGRanges)+100)
maleGRanges<-resize(maleGRanges,fix="start",width=width(maleGRanges)+100)

# make dna string sets
femaleDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, femaleGRanges)
maleDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, maleGRanges)

# get PWMS
hg19_PWMs <- PWMLogn.hg19.MotifDb.Hsap$pwms


################################################################################
# build backgrounds
################################################################################
if(file.exists("/home/og16379/diff_cpg_fm/data/male_background.RData")){
  load("/home/og16379/diff_cpg_fm/data/male_background.RData")
} else{
  male_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=maleDMPseq,
                                              type="logn", algorithm="human",  verbose=TRUE)

  save(male_seq_bg, file="/home/og16379/diff_cpg_fm/data/male_background.RData")
}


if(file.exists("/home/og16379/diff_cpg_fm/data/female_background.RData")){
  load("/home/og16379/diff_cpg_fm/data/female_background.RData")
} else{
  female_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=femaleDMPseq,
                                              type="logn", algorithm="human",verbose=TRUE)

  save(female_seq_bg, file="/home/og16379/diff_cpg_fm/data/female_background.RData")
}



################################################################################
# comparison
################################################################################

# motifs differentially enriched in the maintained sequence (with lognormal background correction)
if(file.exists("/home/og16379/diff_cpg_fm/data/female_male_enrichment.RData")){
  load("/home/og16379/diff_cpg_fm/data/female_male_enrichment.RData")
} else{
  female_male_enrichment <- motifEnrichment(femaleDMPseq, male_seq_bg)
  save(female_male_enrichment, file="/home/og16379/diff_cpg_fm/data/female_male_enrichment.RData")
}

if(file.exists("/home/og16379/diff_cpg_fm/data/male_female_enrichment.RData")){
  load("/home/og16379/diff_cpg_fm/data/male_female_enrichment.RData")
} else{
  male_female_enrichment <- motifEnrichment(maleDMPseq, female_seq_bg)
  save(male_female_enrichment, file="/home/og16379/diff_cpg_fm/data/male_female_enrichment.RData")
}


# group report
# p value set at 0.05
# get unique motifs and remove UW.Motifs
# write to csv
female_male_enrichment_report <- groupReport(female_male_enrichment)
female_male_enrichment_report_pvalue05 <- female_male_enrichment_report[female_male_enrichment_report$p.value < 0.05]
female_male_enrichment_report_pvalue05 <- female_male_enrichment_report_pvalue05[-grep("UW.Motif.",female_male_enrichment_report_pvalue05$target)]
enrichedInDMPsHMInFemales <- unique(female_male_enrichment_report_pvalue05$target)
female_male_enrichment_report_pvalue05df<-as.data.frame(female_male_enrichment_report_pvalue05)
write.csv(female_male_enrichment_report_pvalue05df,"/home/og16379/diff_cpg_fm/data/femaleEnrichedTFBS.csv")

# repeat above for males
male_female_enrichment_report <- groupReport(male_female_enrichment)
male_female_enrichment_report_pvalue05 <- male_female_enrichment_report[male_female_enrichment_report$p.value < 0.05]
male_female_enrichment_report_pvalue05 <- male_female_enrichment_report_pvalue05[-grep("UW.Motif.",male_female_enrichment_report_pvalue05$target)]
enrichedInDMPsHMInMales <- unique(male_female_enrichment_report_pvalue05$target)
male_female_enrichment_report_pvalue05df<-as.data.frame(male_female_enrichment_report_pvalue05)
write.csv(male_female_enrichment_report_pvalue05df,"/home/og16379/diff_cpg_fm/data/maleEnrichedTFBS.csv")

# get list of motifs enriched at sites hypermethylated in females only
femaleOnly<-enrichedInDMPsHMInFemales[-which(enrichedInDMPsHMInFemales %in% enrichedInDMPsHMInMales)]
# get list of motifs enriched at sites hypermethylated in males only
maleOnly<-enrichedInDMPsHMInMales[-which(enrichedInDMPsHMInMales %in% enrichedInDMPsHMInFemales)]
# get list of motifs enriched at sites hypermethylated in females and males
overlap<-enrichedInDMPsHMInMales[which(enrichedInDMPsHMInMales %in% enrichedInDMPsHMInFemales)]


############################################
## complete ##
############################################
