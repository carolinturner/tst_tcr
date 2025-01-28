# TCRseq_Script_2: Calculate diversity indices

library(data.table)
library(tidyverse)
library(ineq)
library(entropy)

# load metadata
meta <- read.csv("data/TCRseq_metadata.csv") %>%
  select(sample,tissue)

# define Simpson function
### unbiased Simpson index, as in E. H. Simpson, Nature 163, 688 (1949); and as
### extended in Tiffeau-Mayer, Phys. Rev. E 109 (2024): https://doi.org/10.1103/PhysRevE.109.064411
## the input is a vector of TCR counts

simpson.sample=function(x){
  n=sum(x,na.rm=TRUE)
  pc=sum((x-1)*x/(n*(n-1)))
  sp_d=1/(pc)
  sp_d
}

# Step 1: full repertoires ####
alpha <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
beta <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7"))

# loop through samples and calculate metrics
for (chain in c("alpha","beta")){
  dat <- if(chain == "alpha"){alpha}else{beta}
  samples <- unique(dat$sample)
  # create empty vectors
  richness <- c()
  ginis <- c()
  shannons.bits <- c()
  simpson <- c()
  sum_expanded <- c()
  sample <- c()
  
  for (i in samples){
    print(paste0("Processing ",i))
    d <- dat %>% filter(sample == i)
    # sum up number of CDR3s independent of gene rearrangement
    d <- d %>% group_by(junction_aa) %>% summarise(total = sum(duplicate_count))
    # rename columns
    colnames(d) <- c("CDR3", "count")
    # calculate metrics
    richness <- c(richness, length(d[[2]])/sum(d[[2]]))
    ginis <- c(ginis, ineq(d[[2]]))
    shannons.bits <- c(shannons.bits, entropy(d[[2]], unit="log2"))
    simpson <- c(simpson, simpson.sample(d[[2]]))
    # subset to expanded clones (count >1)
    exp <- subset(d, count > 1)
    sum_expanded <- c(sum_expanded, sum(exp[[2]]))
    sample <- c(sample,i)
  }
  summary <- data.frame(sample,richness,ginis,shannons.bits,simpson,sum_expanded)
  summary <- merge(summary,meta)
  write.csv(summary,file=paste0("data/Diversity_full-repertoires_",chain,".csv"),row.names = F)
}


# Step 2: down-sampled repertoires ####
alpha <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
beta <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7"))

# loop through samples and calculate metrics
for (chain in c("alpha","beta")){
  dat <- if(chain == "alpha"){alpha}else{beta}
  samples <- unique(dat$sample)
  # create empty vectors
  richness <- c()
  ginis <- c()
  shannons.bits <- c()
  simpson <- c()
  sum_expanded <- c()
  sample <- c()
  
  for (i in samples){
    print(paste0("Processing ",i))
    d <- dat %>% filter(sample == i)
    # sum up number of CDR3s independent of gene rearrangement
    d <- d %>% group_by(junction_aa) %>% summarise(total = sum(duplicate_count))
    # rename columns
    colnames(d) <- c("CDR3", "count")
    # calculate metrics
    richness <- c(richness, length(d[[2]])/sum(d[[2]]))
    ginis <- c(ginis, ineq(d[[2]]))
    shannons.bits <- c(shannons.bits, entropy(d[[2]], unit="log2"))
    simpson <- c(simpson, simpson.sample(d[[2]]))
    # subset to expanded clones (count >1)
    exp <- subset(d, count > 1)
    sum_expanded <- c(sum_expanded, sum(exp[[2]]))
    sample <- c(sample,i)
  }
  summary <- data.frame(sample,richness,ginis,shannons.bits,simpson,sum_expanded)
  summary <- merge(summary,meta)
  write.csv(summary,file=paste0("data/Diversity_down-sampled_",chain,".csv"),row.names = F)
}
