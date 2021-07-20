################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#first, we must normalise our data with the autosomes only
#load packages
suppressMessages(require(bigmelon))
suppressMessages(require(data.table))

#file name
save_file <- 'autosome_norm_beta.xls'
#path to gds file or EPIC file
gds_file <- '/storage/projects/Exeter/UnderSocMeth/USM1/USM_WF.gds'
#open gds file
gfile <- openfn.gds(gds_file)
#get methylated and unmethylated
mns <- methylated(gfile)[,]
uns <- unmethylated(gfile)[,]
#remove sex chromosome probes
autosomes <- !(fData(gfile)$CHR %in% c('X', 'Y'))
print("Dasen normalization on CpGs from autosomes...")
dasen_autosome <- dasen(mns[autosomes, ], uns[autosomes, ], fData(gfile)$Infinium_Design_Type[autosomes])
#write to file
fwrite(data.frame(ID_REF=rownames(dasen_autosome), dasen_autosome, check.names=F), file=save_file, sep='\t')suppressMessages(require(bigmelon))


# next, we can calculate dmps using champ dmp & t test method
# load normalised betas
suppressMessages(require(ChAMP))   ## champ.DMP
suppressMessages(require(data.table))
#arg <- commandArgs(T)
autosome <- TRUE
#autosome <- arg[1]
gds_file <- '/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds'
save_file <- 'diff_cpg_between_F_M.xls'
gfile <- openfn.gds(gds_file)
if(autosome){
    print("Find differential methylated CpG sites on autosomes")
    ## load normalized autosomal beta values
    normed_beta <- '/path/to/file/autosome_norm_beta.xls'
    # betas <- index.gdsn(gfile, 'autobetas')[1:847214,]
    betas <- fread(normed_beta, data.table=FALSE)
    rownames(betas) <- betas$ID_REF
    betas <- betas[, -1]
} else {
    betas <- betas(gfile)[,]
}

## get sex annotation
info <- pData(gfile)[, c('barcode', 'nsex')]
info$nsex <- gsub('2', 'F', info$nsex)
info$nsex <- gsub('1', 'M', info$nsex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1

outlier <- c('200611820013_R08C01', '200603220075_R08C01', '200864580018_R02C01', '200611820020_R07C01')
info <- info[!(info$barcode %in% outlier), ]
betas <- betas[, info$barcode]
pheno_group <- info$nsex


## Find Differential Methylation Positions  by F-test from limma
## add adj.Pval argument and set to 1 if you want all stats for all probes to make the volcano plot in figure 1
f_DMP <- champ.DMP(beta=betas, pheno=pheno_group, arraytype='EPIC')[[1]]

## Find Differential Methylation Positions  by T-test
print("Performing T test: ")
t_Pvals <- apply(betas, 1, function(x){t.test(as.numeric(x)~pheno_group)$p.value})  ## extremly time consuming
t_Padj <- p.adjust(t_Pvals, method = 'BH')
t_DMP <- data.frame(t_Pvals, t_Padj)
t_DMP <- t_DMP[t_DMP$t_Padj < 0.05, ]

#merge the two results together
dmp <- merge(f_DMP, t_DMP, by=0, sort=FALSE)
fwrite(dmp, file=save_file, sep='\t')
