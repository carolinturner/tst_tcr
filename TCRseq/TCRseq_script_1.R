# TCRseq_Script_1: Pre-processing of bulk TCRseq data
# - tidy metadata file
# - tidy HLA data file
# - combine samples in one file
# - down-sample repertoires

# download TCRseq data and metadata from UCL's Research Data Repository: DOI 10.5522/04/28049606
# save all source data in a sub-directory called data

library(tidyverse)
library(data.table)

# Step 1: make tidy metadata and HLA files ####
# meatdata
meta <- read.csv("data/metadata.csv") %>%
  filter(Data_generation_protocol == "UCL_Chain_lab") %>%
  mutate(sample = paste0(UIN,"_",tissue,"_",chain)) %>%
  select(sample,Filename_processed,UIN,chain,tissue,additional_description)
write.csv(meta,"data/TCRseq_metadata.csv",row.names = F)

# hla data
hla <- read.csv("data/TableS3.csv")

# select relevant columns and split HLA gene and allele
dat1 <- hla %>%
  select(UIN, HLA) %>%
  separate_wider_delim(cols = HLA, delim = "_", names = c("HLA","allele")) %>%
  select(-HLA) %>%
  separate_wider_delim(cols = allele, delim = "*", names = c("HLA_gene","allele2"), cols_remove = FALSE) %>%
  select(-allele2) %>%
  separate_wider_delim(cols = allele, delim = ":", names = c("allele", "field2")) %>%
  select(-field2)

# add allele duplicate where second allele not specified
d1 <- dat1 %>%
  group_by(UIN,HLA_gene) %>%
  filter(n() == 1) %>%
  slice(rep(1:n(), each =2)) %>%
  ungroup()
d2 <- dat1 %>%
  group_by(UIN,HLA_gene) %>%
  filter(n() > 1) %>%
  ungroup()
dat2 <- rbind(d1,d2)

# pivot wider
dat3 <- dat2 %>%
  group_by(HLA_gene,UIN) %>%
  mutate(HLA_group = row_number(),
         HLA_id = paste0(HLA_gene,".",HLA_group)) %>%
  ungroup() %>%
  select(-HLA_gene,-HLA_group) %>%
  pivot_wider(id_cols = "UIN", names_from = "HLA_id", values_from = "allele") %>%
  column_to_rownames("UIN")

# combine DQA and DQB alleles (all four combinations)
dat4 <- dat3 %>%
  mutate(DQAB1 = paste0(DQA1.1,"_",DQB1.1),
         DQAB2 = paste0(DQA1.1,"_",DQB1.2),
         DQAB3 = paste0(DQA1.2,"_",DQB1.1),
         DQAB4 = paste0(DQA1.2,"_",DQB1.2)) %>%
  select(-DQA1.1,-DQA1.2,-DQB1.1,-DQB1.2)

# set NA values (where no HLA allele imputation is provided) to MISSING
dat5 <- replace(dat4, is.na(dat4),"MISSING")

# save to file
write.csv(dat5,"data/hladata.csv",row.names = T)

# Step 2: combine full repertoires in one file per chain ####

# split metadata by chain
meta.a <- meta %>% filter(chain == "alpha") 
meta.b <- meta %>% filter(chain == "beta")

# combine repertoires and add relevant metadata annotations
for (chain in c("alpha","beta")){
  m <- if(chain == "alpha"){meta.a}else{meta.b}
  summary <- data.frame()
  for (i in 1:nrow(m)){
    d <- fread(paste0("data/",m[i,"Filename_processed"]))
    d <- subset(d, productive == "T") # keep only productive rearrangements
    d$UIN <- m[i,"UIN"]
    d$tissue <- m[i,"tissue"]
    d$sample <- m[i,"sample"]
    d$replicate <- m[i,"additional_description"]
    d$bioidentity <- paste0(d$v_call,d$junction_aa,d$j_call)
    summary <- rbind(summary,d)
  }
  fwrite(summary,file = paste0("data/combined_",chain,".csv.gz"))
}

# Step 3: down-sampling of blood and TST samples ####

# load relevant repertoires
alpha <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
beta <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7"))

# make tables of repertoire sizes to manually look for suitable sub-sample size
repsizes.a <- alpha %>%
  group_by(tissue, sample) %>%
  summarise(total = sum(duplicate_count))
repsizes.b <- beta %>%
  group_by(tissue, sample) %>%
  summarise(total = sum(duplicate_count))

# define sub-sample size
num <- 16000
spl.a <- repsizes.a %>% filter(total >= num)
spl.b <- repsizes.b %>% filter(total >= num)

# define output path
dir.create(paste0("data/Subsampling_",num))
path_out <- paste0("data/Subsampling_",num,"/")

# extract list of blood and skin samples that meet down-sampling size
meta.ss <- meta %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7") & (sample %in% spl.a$sample | sample %in% spl.b$sample))

# sub-sample each file
for (i in 1:nrow(meta.ss)){
  d <- fread(paste0("data/",meta.ss[i,"Filename_processed"]))
  d <- subset(d, productive == "T") # keep only productive rearrangements
  d <- d %>% uncount(duplicate_count) # un-count
  d$count <- 1 # add count column
  # take a random sample
  set.seed(123)
  subset <- sample(1:nrow(d),num)
  sample <- d[subset,]
  # sum up original sequence ids
  counts <- sample %>% group_by(sequence_id) %>% summarise(duplicate_count=sum(count))
  sample$count <- NULL
  # add all columns
  setDT(sample)
  setDT(counts)
  full <- sample[counts, mult="first",on="sequence_id",nomatch=0L]
  # write to file
  write.csv(full,file=paste0(path_out,"subsampled_",meta.ss[i,"sample"],".csv"),row.names = F)
}

# Step 4: combine down-sampled repertoires in one file per chain ####

# split metadata by chain
meta.a <- meta.ss %>% filter(chain == "alpha")
meta.b <- meta.ss %>% filter(chain == "beta")

# combine repertoires
for (chain in c("alpha","beta")){
  m <- if(chain == "alpha"){meta.a}else{meta.b}
  summary <- data.frame()
  for (i in 1:nrow(m)){
    d <- read.csv(paste0(path_out,"subsampled_",m[i,"sample"],".csv"))
    d$UIN <- m[i,"UIN"]
    d$tissue <- m[i,"tissue"]
    d$sample <- m[i,"sample"]
    d$bioidentity <- paste0(d$v_call,d$junction_aa,d$j_call)
    summary <- rbind(summary,d)
  }
  fwrite(summary,file = paste0("data/combined_subsampled_",chain,".csv.gz"))
}
