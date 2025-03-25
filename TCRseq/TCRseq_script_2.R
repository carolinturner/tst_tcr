# TCRseq_Script_2: Calculate diversity indices

library(data.table)
library(tidyverse)
library(ineq)
library(entropy)

# load metadata
meta <- read.csv("data/TCRseq_metadata.csv") %>%
  select(sample,tissue)

# define inverse Simpson function
### population Simpson index, as in E. H. Simpson, Nature 163, 688 (1949)
## the input is a vector of TCR counts

invsimpson.population=function(x){
  n=sum(x,na.rm=TRUE)
  p=x/n
  pc=sum(p^2)
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
  gini <- c()
  shannon <- c()
  invsimpson <- c()
  sum_expanded <- c()
  sample <- c()
  
  for (i in samples){
    print(paste0("Processing ",i))
    d <- dat %>% filter(sample == i)
    richness <- c(richness, length(d$duplicate_count))
    gini <- c(gini, ineq(d$duplicate_count))
    shannon <- c(shannon, exp(entropy(d$duplicate_count, method="ML")))
    invsimpson <- c(invsimpson, invsimpson.population(d$duplicate_count))
    
    # subset to expanded clones (count >1)
    exp <- subset(d, duplicate_count > 1)
    sum_expanded <- c(sum_expanded, sum(exp$duplicate_count))
    sample <- c(sample,i)
  }
  summary <- data.frame(sample,richness,gini,shannon,invsimpson,sum_expanded)
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
  gini <- c()
  shannon <- c()
  invsimpson <- c()
  sum_expanded <- c()
  sample <- c()
  
  for (i in samples){
    print(paste0("Processing ",i))
    d <- dat %>% filter(sample == i)
    richness <- c(richness, length(d$duplicate_count))
    gini <- c(gini, ineq(d$duplicate_count))
    shannon <- c(shannon, exp(entropy(d$duplicate_count, method="ML")))
    invsimpson <- c(invsimpson, invsimpson.population(d$duplicate_count))
    
    # subset to expanded clones (count >1)
    exp <- subset(d, duplicate_count > 1)
    sum_expanded <- c(sum_expanded, sum(exp$duplicate_count))
    sample <- c(sample,i)
  }
  summary <- data.frame(sample,richness,gini,shannon,invsimpson,sum_expanded)
  summary <- merge(summary,meta)
  write.csv(summary,file=paste0("data/Diversity_down-sampled_",chain,".csv"),row.names = F)
}
