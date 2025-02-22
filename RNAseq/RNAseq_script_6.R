library(tidyverse)
library(XGR)
library(GeneOverlap)

# load annotations file to allow mapping of Ensembl Gene IDs to Gene Symbols
annotations <- read.csv("annotations_ref111.csv", header=T)
annotations <- annotations %>% 
  select(ensembl_gene_id, external_gene_name) %>%
  dplyr::rename(Id = ensembl_gene_id, GeneSymbol = external_gene_name)

# load genes of interest (= output from DeSeq2 analysis)
dat <- read.table("SARtools_output_D7vsD2/D7vsD2.complete.txt", header = T)
dat1 <- dat %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  select(Id) %>%
  left_join(annotations) %>%
  select(-Id) %>%
  distinct(GeneSymbol) %>%
  mutate(Group = "Day 2 TST")
dat2 <- dat %>%
  filter(log2FoldChange >=1 & padj < 0.05) %>%
  select(Id) %>%
  left_join(annotations) %>%
  select(-Id) %>%
  distinct(GeneSymbol) %>%
  mutate(Group = "Day 7 TST")

gene_list <- rbind(dat1,dat2)

# load background genes (= all genes that were used as input for differential gene expression)
background <- dat %>%
  select(Id) %>%
  left_join(annotations)
  
background_genes <- unique(background$GeneSymbol)

# Pathways analysis in XGR
ontology <- "MsigdbC2REACTOME" ## choose ontology, here Reactome

pathways <- gene_list %>% 
  group_by(Group) %>% 
  do(as.data.frame(tryCatch(filter(xEnrichViewer(xEnricherGenes(data=.$GeneSymbol,ontology=ontology,background=background_genes,p.adjust.method="BH",check.symbol.identity
                                                                =TRUE), sortBy = "fdr", top_num=100, details=T), adjp < 0.05),error=function(e)NA)))

# Clustering:
mydf <- pathways %>% 
  select(name,members_Overlap)

make_list <- function(variable_genes){
  as.list(strsplit(variable_genes,", "))
}
mydf_d <- apply(mydf,1,make_list)
names(mydf_d) <- mydf$name

mydf_list_unique <- lapply(mydf_d, unique)

mylist <- rep(list(),length(mydf_d))

for(i in 1:length(mydf_d)) {
  
  mylist[[i]] <- (mydf_d[[i]]$members_Overlap) ## overlap vs annotation
}

names(mylist) <- mydf$name

# Calculate pairwise Jaccard index between any two pathways
result <- sapply(mylist,function(x) 
  sapply(mylist,function(y,x)length(intersect(x,y))/length(union(x,y)),x))

# Hierarchical clustering of Jaccard indices
clusters <- hclust(dist(result))


clusterCut <- cutree(clusters, 20) ## change here the desired number of clusters (here 20) 
clusterCut <- as.data.frame(clusterCut)
clusterCut$group <- clusterCut$clusterCut
clusterCut$name <- row.names(clusterCut)
new_list_grouped <- left_join(pathways,clusterCut)

# write to file
write.csv(new_list_grouped, "pathway_analysis.csv",row.names = F)
