## RNAseq_script_1: Prepare data and differential gene expression analysis

# download data from ArrayExpress (https://www.ebi.ac.uk/arrayexpress), accession number E-MTAB-14687 
# metadata = sdrf.tsv
# counts data = rawcounts.csv
# tpm data = tpm_PC0.001_log2_genesymbol_dedup.csv

library(tidyverse)

# make tidy metadata file
meta <- read.table("data/sdrf.tsv",header=TRUE,fill = TRUE,sep = "\t")
meta <- meta %>%
  select(Source.Name,Characteristics.individual.,Characteristics.stimulus.) %>%
  rename(sample = Source.Name,
         UIN = Characteristics.individual.,
         Stimulant = Characteristics.stimulus.)
write.csv(meta,"data/RNAseq_metadata.csv",row.names = F)

# make tidy annotation file
data <- read.csv("data/rawcounts.csv")
annot <- data %>% select(ensembl_gene_id,external_gene_name,gene_biotype)
write.csv(annot,"data/annotations_ref111.csv",row.names = F)

# filter rawcounts matrix to exclude genes annotated as any type of pseudogene
biotype_all <- unique(annot$gene_biotype)
biotype_clean <- biotype_all[-grep("*pseudogene",biotype_all)]

d <- data %>% filter(gene_biotype %in% biotype_clean) %>% select(-gene_biotype,-external_gene_name)
write.csv(d,"data/rawcounts_biotype_clean.csv",row.names = F)

# split rawcounts matrix as input for SARtools
d <- read.csv("data/rawcounts_biotype_clean.csv")
gene_ids <- as.data.frame(d[,1])
d1 <- round(d[,2:ncol(d)]) # round count values to integers
names <- colnames(d1)

# make a new folder where single count files will be saved to in .txt format
dir.create("data/raw")

# loop through samples, merge with list of genes and save as .txt file with sample name
for (i in 1:length(names)){
  sample <- d1[,i]
  merged <- cbind(gene_ids,sample)
  write.table(merged,file=paste0("data/raw/",names[i],".txt"),row.names=FALSE,col.names=FALSE,sep="\t")
}

# make metadata tables for each comparison of interest
meta <- read.csv("data/RNAseq_metadata.csv")
m1 <- meta %>%
  filter(Stimulant %in% c("saline","TST_D2")) %>%
  mutate(Stimulant = recode(Stimulant,
                            "TST_D2" = "TSTd2"),
         Filename = paste0(sample,".txt")) %>%
  select(sample,Filename,Stimulant)
write.table(m1,"data/SARtools_meta_TSTd2vsSal.txt",sep="\t",quote = F, row.names = F)

m2 <- meta %>%
  filter(Stimulant %in% c("saline","TST_D7")) %>%
  mutate(Stimulant = recode(Stimulant,
                            "TST_D7" = "TSTd7"),
         Filename = paste0(sample,".txt")) %>%
  select(sample,Filename,Stimulant)
write.table(m2,"data/SARtools_meta_TSTd7vsSal.txt",sep="\t",quote = F, row.names = F)

m3 <- meta %>%
  mutate(Stimulant = recode(Stimulant,
                            "TST_D7" = "TST",
                            "TST_D2" = "TST"),
         Filename = paste0(sample,".txt")) %>%
  select(sample,Filename,Stimulant)
write.table(m3,"data/SARtools_meta_TSTvsSal.txt",sep="\t",quote=F, row.names = F)

m4 <- meta %>%
  filter(Stimulant %in% c("TST_D2","TST_D7")) %>%
  mutate(Stimulant = recode(Stimulant,
                            "TST_D2" = "D2",
                            "TST_D7" = "D7"),
         Filename = paste0(sample,".txt")) %>%
  select(sample,Filename,Stimulant)
write.table(m4,"data/SARtools_meta_D7vsD2.txt",sep="\t",quote = F, row.names = F)
