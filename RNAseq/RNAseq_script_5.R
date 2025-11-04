## RNAseq_script_5: SARtools DeSeq2 analysis D7 vs D2 (within integrated TST response)

# split rawcounts matrix as input for SARtools
d <- read.csv("data/rawcounts_integrated-TST-transcriptome.csv")
gene_ids <- as.data.frame(d[,1])
d1 <- round(d[,2:ncol(d)]) # round count values to integers
names <- colnames(d1)

# make a new folder where single count files will be saved to in .txt format
dir.create("data/raw_integrated")

# loop through samples, merge with list of genes and save as .txt file with sample name
for (i in 1:length(names)){
  sample <- d1[,i]
  merged <- cbind(gene_ids,sample)
  write.table(merged,file=paste0("data/raw_integrated/",names[i],".txt"),row.names=FALSE,col.names=FALSE,sep="\t")
}


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

workDir <- "C:/path/to/your/working/directory/data/" # working directory for the R session

projectName <- "D7 vs D2"                            # name of the project
author <- "Your name"                                # author of the statistical analysis/report

targetFile <- "SARtools_meta_D7vsD2.txt"             # path to the design/target file
rawDir <- "raw_integrated"                           # path to the directory containing raw counts files
featuresToRemove <- NULL                             # names of the features to be removed
# (specific HTSeq-count information and rRNA for example)
# NULL if no feature to remove

varInt <- "Stimulant"                                # factor of interest
condRef <- "D2"                                      # reference biological condition
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
dir.create("SARtools_output_D7vsD2")
file.rename("./tables/D7vsD2.complete.txt","SARtools_output_D7vsD2/D7vsD2.complete.txt")
file.rename("./tables/D7vsD2.up.txt","SARtools_output_D7vsD2/D7vsD2.up.txt")