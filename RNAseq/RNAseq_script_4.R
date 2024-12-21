## RNAseq_script_4: SARtools DeSeq2 analysis all TST vs Saline

# SARtools template script downloaded from:
# https://github.com/PF2-pasteur-fr/SARTools/blob/master/template_script_DESeq2.r

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### March 20th, 2018
### designed to be executed with SARTools 1.6.7
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "C:/path/to/your/working/directory/"      # working directory for the R session

projectName <- "TST vs Saline"                       # name of the project
author <- "Your name"                                # author of the statistical analysis/report

targetFile <- "SARtools_meta_TSTvsSal.txt"           # path to the design/target file
rawDir <- "raw"                                      # path to the directory containing raw counts files
featuresToRemove <- NULL                             # names of the features to be removed
# (specific HTSeq-count information and rRNA for example)
# NULL if no feature to remove

varInt <- "Stimulant"                                # factor of interest
condRef <- "saline"                                  # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- FALSE                                 # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "fdr"                               # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

forceCairoGraph <- TRUE

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
#install.packages("devtools")
#devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
#majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
#exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
#save.image(file=paste0(projectName, ".RData"))

# generating HTML report
#writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
#                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
#                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
#                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
#                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
#                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# move relevant output tables to avoid overwriting by additional SARtools analyses
dir.create("SARtools_output_TST")
file.rename("./tables/TSTvssaline.complete.txt","SARtools_output_TST/TSTvssaline.complete.txt")
file.rename("./tables/TSTvssaline.up.txt","SARtools_output_TST/TSTvssaline.up.txt")

# up-regulated genes (FDR<0.05, log2FC>=1)
up <- read.table("SARtools_output_TST/TSTvssaline.up.txt", header=T) %>% 
  filter(log2FoldChange >= 1) %>%
  select(Id)

# subset tpm matrix to keep only genes up-regulated in TST
tpm <- read.csv("tpm_PC0.001_log2_genesymbol_dedup.csv")
annot <- read.csv("annotations_ref111.csv",na.strings = c("",NA)) %>%
  select(-gene_biotype)
tpm.up <- up %>%
  dplyr::rename("ensembl_gene_id" = "Id") %>%
  left_join(annot) %>%
  select(external_gene_name) %>%
  distinct() %>%
  na.omit() %>%
  left_join(tpm)
write.csv(tpm.up, file="tpm_TST-transcriptome.csv",row.names = F)

# make input file for IPA upstream regulator analysis
up <- read.table("SARtools_output_TST/TSTvssaline.up.txt", header=T) %>% 
  filter(log2FoldChange >= 1) %>%
  select(Id,log2FoldChange,padj)
write.table(up,"IPA_input_allTST-vs-sal.txt",sep="\t",quote=F,row.names = F)