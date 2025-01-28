# TCRseq_script_7: Calculate abundance of private, antigen-reactive TCRs in blood and TST samples

library(data.table)
library(tidyverse)

# load metadata and identify UINs with in-vitro stimulation data
meta <- read.csv("data/TCRseq_metadata.csv")
id <- meta %>% filter(tissue == "PBMC_PPD") %>% pull(UIN) %>% unique()

# Step 1: full repertoires ####
alpha <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7") & UIN %in% id)
beta <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7") & UIN %in% id)

# loop through different expansion thresholds and quantify published CDR3s
for (f in c(0:4)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs
  dat.a <- alpha %>%
    group_by(sample,junction_aa) %>%
    mutate(cdr3_count = sum(duplicate_count)) %>%
    filter(cdr3_count >f) %>%
    ungroup()
  dat.b <- beta %>%
    group_by(sample,junction_aa) %>%
    mutate(cdr3_count = sum(duplicate_count)) %>%
    filter(cdr3_count >f) %>%
    ungroup()
  
  # make empty results data frame
  results <- data.frame()
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    dat <- if(chain == "alpha"){dat.a}else{dat.b}
    # load references
    PPD.CDR3 <- as.character(read.csv(paste0("data/PPD_",chain,".csv"), header = F)$V1)
    TT.CDR3 <- as.character(read.csv(paste0("data/TT_",chain,".csv"), header = F)$V1)
    
    # annotate data with Ag-reactivity
    dat$PPD <- ifelse(dat$junction_aa %in% PPD.CDR3, "PPD", NA)
    dat$TT <- ifelse(dat$junction_aa %in% TT.CDR3, "TT", NA)
    
    # summary dataframe
    summary <- dat %>% filter(!if_all(c(PPD,TT), is.na))

    # sum total and unique counts for each sample
    counts <- dat %>%
      group_by(tissue,sample) %>%
      summarise(total = sum(duplicate_count),
                unique = length(duplicate_count))
    
    # sum Ag-reactive counts for each sample
    PPD.counts <- summary %>%
      filter(PPD == "PPD") %>%
      group_by(tissue,sample) %>%
      summarise(ag.total = sum(duplicate_count),
                ag.unique = length(duplicate_count))
    TT.counts <- summary %>%
      filter(TT == "TT") %>%
      group_by(tissue,sample) %>%
      summarise(ag.total = sum(duplicate_count),
                ag.unique = length(duplicate_count))
    
    # merge data frames and calculate percentage
    PPD.dat <- left_join(counts,PPD.counts) %>%
      mutate(pct = ag.total/total*100,
             pct.unique = ag.unique/unique*100,
             pct = replace_na(pct, 0),
             pct.unique = replace_na(pct.unique, 0),
             Antigen = "PPD") %>%
      select(tissue,sample,Antigen,pct,pct.unique)
    TT.dat <- left_join(counts,TT.counts) %>%
      mutate(pct = ag.total/total*100,
             pct.unique = ag.unique/unique*100,
             pct = replace_na(pct, 0),
             pct.unique = replace_na(pct.unique, 0),
             Antigen = "TT") %>%
      select(tissue,sample,Antigen,pct,pct.unique)
    merge <- rbind(PPD.dat,TT.dat)
    merge$chain <- chain
    results <- rbind(results,merge)
  }
  write.csv(results,paste0("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr",f,".csv"),row.names = F)
}

# Step 2: down-sampled repertoires ####
alpha <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7") & UIN %in% id)
beta <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7") & UIN %in% id)

# loop through different expansion thresholds and quantify published CDR3s
for (f in c(0:4)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs
  dat.a <- alpha %>%
    group_by(sample,junction_aa) %>%
    mutate(cdr3_count = sum(duplicate_count)) %>%
    filter(cdr3_count >f) %>%
    ungroup()
  dat.b <- beta %>%
    group_by(sample,junction_aa) %>%
    mutate(cdr3_count = sum(duplicate_count)) %>%
    filter(cdr3_count >f) %>%
    ungroup()
  
  # make empty results data frame
  results <- data.frame()
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    dat <- if(chain == "alpha"){dat.a}else{dat.b}
    # load references
    PPD.CDR3 <- as.character(read.csv(paste0("data/PPD_",chain,".csv"), header = F)$V1)
    TT.CDR3 <- as.character(read.csv(paste0("data/TT_",chain,".csv"), header = F)$V1)
    
    # annotate data with Ag-reactivity
    dat$PPD <- ifelse(dat$junction_aa %in% PPD.CDR3, "PPD", NA)
    dat$TT <- ifelse(dat$junction_aa %in% TT.CDR3, "TT", NA)
    
    # summary dataframe
    summary <- dat %>% filter(!if_all(c(PPD,TT), is.na))

    # sum total and unique counts for each sample
    counts <- dat %>%
      group_by(tissue,sample) %>%
      summarise(total = sum(duplicate_count),
                unique = length(duplicate_count))
    
    # sum Ag-reactive counts for each sample
    PPD.counts <- summary %>%
      filter(PPD == "PPD") %>%
      group_by(tissue,sample) %>%
      summarise(ag.total = sum(duplicate_count),
                ag.unique = length(duplicate_count))
    TT.counts <- summary %>%
      filter(TT == "TT") %>%
      group_by(tissue,sample) %>%
      summarise(ag.total = sum(duplicate_count),
                ag.unique = length(duplicate_count))
    
    # merge data frames and calculate percentage
    PPD.dat <- left_join(counts,PPD.counts) %>%
      mutate(pct = ag.total/total*100,
             pct.unique = ag.unique/unique*100,
             pct = replace_na(pct, 0),
             pct.unique = replace_na(pct.unique, 0),
             Antigen = "PPD") %>%
      select(tissue,sample,Antigen,pct,pct.unique)
    TT.dat <- left_join(counts,TT.counts) %>%
      mutate(pct = ag.total/total*100,
             pct.unique = ag.unique/unique*100,
             pct = replace_na(pct, 0),
             pct.unique = replace_na(pct.unique, 0),
             Antigen = "TT") %>%
      select(tissue,sample,Antigen,pct,pct.unique)
    merge <- rbind(PPD.dat,TT.dat)
    merge$chain <- chain
    results <- rbind(results,merge)
  }
  write.csv(results,paste0("data/Summary_down-sampled_private-Ag-abundance_expanded_gr",f,".csv"),row.names = F)
}

