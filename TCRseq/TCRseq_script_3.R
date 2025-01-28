# TCRseq_Script_3: Calculate abundance of published antigen-reactive TCRs in blood and TST samples

library(data.table)
library(tidyverse)

# load published sequences
ref <- read.csv("data/TableS2.csv")

# Step 1: full repertoires ####
alpha <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
beta <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7"))

# loop through different expansion thresholds and quantify published CDR3s
for (f in c(0:1)) {
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
  
  for (Chain in c("alpha","beta")){
    print(paste0("Chain: ",Chain))
    # assign dataset
    dat <- if(Chain == "alpha"){dat.a}else{dat.b}
    # load references
    CMV <- ref %>% filter(reactivity == "CMV" & chain == Chain) %>% pull(CDR3) %>% unique()
    EBV <- ref %>% filter(reactivity == "EBV" & chain == Chain) %>% pull(CDR3) %>% unique()
    Mtb <- ref %>% filter(reactivity == "Mtb" & chain == Chain) %>% pull(CDR3) %>% unique()
    
    # annotate data with Ag-reactivity
    dat$CMV <- ifelse(dat$junction_aa %in% CMV, "CMV", NA)
    dat$EBV <- ifelse(dat$junction_aa %in% EBV, "EBV", NA)
    dat$Mtb <- ifelse(dat$junction_aa %in% Mtb, "Mtb", NA)
    
    # summary dataframe
    summary <- dat %>% filter(!if_all(c(CMV,EBV,Mtb), is.na))
    
    # sum total counts for each sample
    dat <- dat %>%
      group_by(tissue,sample) %>%
      summarise(count = sum(duplicate_count))
    
    # sum Ag-reactive counts for each sample
    CMV.dat <- summary %>%
      filter(CMV == "CMV") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    EBV.dat <- summary %>%
      filter(EBV == "EBV") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    Mtb.dat <- summary %>%
      filter(Mtb == "Mtb") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    
    # merge dataframes and calculate percentage
    merge.CMV <- left_join(dat,CMV.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "CMV") %>%
      select(tissue,Antigen,pct)
    merge.EBV <- left_join(dat,EBV.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "EBV") %>%
      select(tissue,Antigen,pct)
    merge.Mtb <- left_join(dat,Mtb.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "Mtb") %>%
      select(tissue,Antigen,pct)
    
    # save summary file
    summary <- rbind(merge.CMV,merge.EBV,merge.Mtb)
    write.csv(summary,paste0("data/Published-Ag-abundance_full-repertoires_expanded_gr",f,"_",Chain,".csv"),row.names = F)
  }
}

# Step2: down-sampled repertoires ####

alpha <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
beta <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood", "TST_D2", "TST_D7"))

# loop through different expansion thresholds and quantify published CDR3s
for (f in c(0:1)) {
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
  
  for (Chain in c("alpha","beta")){
    print(paste0("Chain: ",Chain))
    # assign dataset
    dat <- if(Chain == "alpha"){dat.a}else{dat.b}
    # load references
    CMV <- ref %>% filter(reactivity == "CMV" & chain == Chain) %>% pull(CDR3) %>% unique()
    EBV <- ref %>% filter(reactivity == "EBV" & chain == Chain) %>% pull(CDR3) %>% unique()
    Mtb <- ref %>% filter(reactivity == "Mtb" & chain == Chain) %>% pull(CDR3) %>% unique()
    
    # annotate data with Ag-reactivity
    dat$CMV <- ifelse(dat$junction_aa %in% CMV, "CMV", NA)
    dat$EBV <- ifelse(dat$junction_aa %in% EBV, "EBV", NA)
    dat$Mtb <- ifelse(dat$junction_aa %in% Mtb, "Mtb", NA)
    
    # summary dataframe
    summary <- dat %>% filter(!if_all(c(CMV,EBV,Mtb), is.na))

    # sum total counts for each sample
    dat <- dat %>%
      group_by(tissue,sample) %>%
      summarise(count = sum(duplicate_count))
    
    # sum Ag-reactive counts for each sample
    CMV.dat <- summary %>%
      filter(CMV == "CMV") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    EBV.dat <- summary %>%
      filter(EBV == "EBV") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    Mtb.dat <- summary %>%
      filter(Mtb == "Mtb") %>%
      group_by(tissue,sample) %>%
      summarise(ag.count = sum(duplicate_count))
    
    # merge dataframes and calculate percentage
    merge.CMV <- left_join(dat,CMV.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "CMV") %>%
      select(tissue,Antigen,pct)
    merge.EBV <- left_join(dat,EBV.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "EBV") %>%
      select(tissue,Antigen,pct)
    merge.Mtb <- left_join(dat,Mtb.dat) %>%
      mutate(pct = ag.count/count*100,
             pct = replace_na(pct,0),
             Antigen = "Mtb") %>%
      select(tissue,Antigen,pct)
    
    # save summary file
    summary <- rbind(merge.CMV,merge.EBV,merge.Mtb)
    write.csv(summary,paste0("data/Published-Ag-abundance_down-sampled_expanded_gr",f,"_",Chain,".csv"),row.names = F)
  }
}
