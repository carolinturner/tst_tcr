## RNAseq_script_8: Identification of statistically significant co-regulated upstream regulator modules (TST_D7 vs saline)

# script adapted from Blanca Sanz-Magallon Duque De Estrada

library(tidyverse)

## A. Input files - replace as necessary
tpm <- read.csv("data/tpm_D7-transcriptome.csv", row.names = 1)
ipa <- read.table("data/IPA_output/TSTd7vsSaline.txt", sep = "\t", skip=2, header=T)
meta <- read.csv("data/RNAseq_metadata.csv", header = TRUE, row.names = 1)

## B. Create a co-correlation matrix of the TST transcriptome (including only relevant TST samples)
meta.ss <- meta %>% filter(Stimulant %in% c("TST_D7"))
tpm.ss <- tpm[,rownames(meta.ss)]

TST_transcriptome <- tpm.ss
TST_transcriptome <- TST_transcriptome[order(rownames(TST_transcriptome)),]
r_tst <- cor(t(TST_transcriptome), method = "spearman")

## C. Format the IPA output file
types <- c("cytokine", "transcription regulator", "kinase", "transmembrane receptor")
ipa <- ipa %>%
  filter(B.H.corrected.p.value < 0.05 & Predicted.Activation.State == "Activated" & Molecule.Type %in% types) %>%
  select(Upstream.Regulator, Molecule.Type, B.H.corrected.p.value, Target.Molecules.in.Dataset)
colnames(ipa) <- c("regulator","type","FDR","target")
ipa <- tibble::rowid_to_column(ipa, "ID") #add an ID column corresponding to each cluster
ipa$target <- as.character(ipa$target) #target column is integer -> convert to character for next step
ipa_sep <- separate_rows(ipa, target, sep = ",|/", convert = TRUE) #separate clusters into individual rows
ipa_sep_reduced <- subset(ipa_sep, target %in% colnames(r_tst)) #only the genes found in the TST transcriptome co-correlation matrix
ipa_sep_reduced %>%
  group_by(ID) %>%
  mutate("count" = n()) -> ipa_sep_reduced # add a column with the cluster size

## D. Calculate the average of all pairwise correlations for the target genes of each upstream regulator cluster
avg_interactomes <- matrix(nrow=0, ncol=2) # empty data frame to add the calculated average correlation values
colnames(avg_interactomes) <- c("regulator", "avg_corr")
for (i in 1:max(ipa_sep_reduced$ID)){
  ipa_sep_reduced %>%
    filter(ID == i) -> hub
  genes <- hub$target[which(hub$target %in% row.names(r_tst))]
  r_tst_network <- r_tst[genes,genes]
  lower_matrix <- lower.tri(r_tst_network)
  r_tst_network_lower <- r_tst_network[lower_matrix]
  mean_interactome <- mean(r_tst_network_lower[r_tst_network_lower>=0])
  avg_interactomes <- rbind(avg_interactomes, c(as.character(unique(hub$regulator)), mean_interactome))
}
# add the calculated cluster average correlations to a data frame containing the other information
ipa_sep_reduced$target <- NULL
ipa_sep_reduced <- unique(ipa_sep_reduced)
ipa_sep_reduced <- merge(ipa_sep_reduced, avg_interactomes, by.x = "regulator", by.y = "regulator") # merge the ipa output file and the correlation values by common regulator name

## E. Generate frequency distributions of the expression of randomly selected groups of genes of a range of sizes
## notes:
##   p = the maximum size of the group of genes to test (ie. range 4 - k). 4 is the lowest number for the range, the correlation function doesn't work for size < 4
##   r = the number of random samples used to calculate the frequency distribution for a given size -- the larger the number the more accurate, but at the cost of longer running time
p <- 578 # check ipa_sep_reduced for largest count
r <- 100
gene_names <- rownames(TST_transcriptome)

random_distributions_100 <- matrix(nrow = 0, ncol = 6) # empty matrix to fill with the distribution values
colnames(random_distributions_100) <- c("size", "mean", "sd", "84.13%", "97.72%",	"99.87%")
for (k in (4:p)){
  random_genes <- matrix(nrow = 0, ncol = 2)
  colnames(random_genes) <- c("iteration", "average_correlation")
  for (i in (1:r)){ # Do this 100 times (could also try doing 1000 times if it doesn't take too long for the loop to run)
    genes <- sample(gene_names, k) # random sample from TST transcriptome (note: no seed set, so repeat runs will yield slightly different results)
    r_tst_network <- r_tst[genes,genes]
    lower_matrix <- lower.tri(r_tst_network)
    r_tst_network_lower <- r_tst_network[lower_matrix]
    mean_interactome <- mean(r_tst_network_lower[r_tst_network_lower>=0]) # Calculate the average correlation for each of the 100 (or 1000) random clusters
    random_genes <- rbind(random_genes, c(i, mean_interactome))
  }
  mean <- mean(as.numeric(random_genes[,2]))
  sd <- sd(as.numeric(random_genes[,2]))
  random_distributions_100 <- rbind(random_distributions_100, c(k, mean, sd, mean + sd, mean + 2*sd, mean + 3*sd))
}
colnames(random_distributions_100) <- c("size", "mean", "sd", "X84.13", "X97.72",	"X99.86")
random_distributions_100_df <- as.data.frame(random_distributions_100)

## F. Determining which clusters are FDR significant
ipa_avg_corr <- ipa_sep_reduced[ipa_sep_reduced$count >= 4,] # take out clusters smaller than 4
ipa_avg_corr$rc_mean <- random_distributions_100_df$mean[match(ipa_avg_corr$count, random_distributions_100_df$size)] #mean
ipa_avg_corr$rc_sd <- random_distributions_100_df$sd[match(ipa_avg_corr$count, random_distributions_100_df$size)] #sd
ipa_avg_corr$avg_corr <- as.numeric(ipa_avg_corr$avg_corr)
ipa_avg_corr$zscore <- (ipa_avg_corr$avg_corr - ipa_avg_corr$rc_mean)/ipa_avg_corr$rc_sd #zscores
ipa_avg_corr$pvalue <- pnorm(ipa_avg_corr$zscore, lower.tail = FALSE) #p-values
ipa_less_0.05 <- ipa_avg_corr[which(ipa_avg_corr$pvalue <= 0.05),] # p-values <= 0.05
fdr_threshold <- nrow(ipa_less_0.05)
ipa_less_0.05$adj_pvalue <- ipa_less_0.05$pvalue*fdr_threshold # calculating the adjusted p-value
ipa_significant <- ipa_less_0.05[which(ipa_less_0.05$adj_pvalue <= 0.05),] # selecting clusters with adj. p-value <= 0.05

## G. Output file
write.csv(ipa_significant, "data/URA_TSTd7vsSaline_fdrsig_578size.csv", row.names = F) # clusters with adj. p-value <= 0.05

# prepare input files for Gephi
# add correlation analysis outcome to filtered IPA results
ipa_significant <- ipa_significant %>% 
  select(regulator, adj_pvalue)
ipa <- ipa %>% left_join(ipa_significant, by = "regulator")
ipa$significant.in.correl.analysis <- ifelse(ipa$adj_pvalue < 0.05, "Yes", "No")
ipa <- ipa %>% select(-adj_pvalue)  

# select only significant upstream regulators
dat <- subset(ipa, significant.in.correl.analysis == "Yes")
dat <- dat %>% select(regulator,type,target,FDR)

# transform p-values so that most significant regulator is associated with largest value
dat$FDR <- log10(dat$FDR) * (-1)

########################
## Make a 'Node' file ##
########################

# make an upstream regulator and gene vector
UR <- dat$regulator
gene <- unique(unlist(strsplit(dat$target, split = ",")))

# make weight vectors
weight.UR <- dat$FDR
weight.gene <- rep(1, times = length(gene))

# make class vectors for upstream regulators and genes
class.UR <- dat$type
class.gene <- rep("Gene", times = length(gene)) 

# make a 'node' data frame and save to .csv file
node.df <- data.frame(Id = c(UR,gene),
                      Label = c(UR, gene),
                      Class = c(class.UR, class.gene),
                      Weight = c(weight.UR, weight.gene))

write.csv(node.df, "data/Nodes_TSTd7vsSaline.csv", row.names = F)

#########################
## Make an 'Edge' file ##
#########################

# drop FDR and type column
dat <- dat %>% select(-c(FDR, type))

# make a separate row for each combination of upstream regulator and target molecule
edge <- tidyr::separate_rows(dat, target, sep=",")

# add a column specifying the edge type
edge$Type <- "undirected"

# change column names and save to file
colnames(edge) <- c("Target", "Source", "Type")

write.csv(edge, "data/Edges_TSTd7vsSaline.csv", row.names = F)
