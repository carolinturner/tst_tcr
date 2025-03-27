library(tidyverse)
library(data.table)

# load metadata and select samples with in vitro stim data
meta <- read.csv("data/TCRseq_metadata.csv") 
id <- meta %>% filter(tissue == "PBMC_PPD") %>% pull(UIN) %>% unique()

### Step 1: private antigen-reactive CDR3s per individual (beta only) ####
for (i in id){
  nil <- read.csv(paste0("data/Summary_Nil-reactive-CDR3s_exp-thr3_beta.csv")) %>%
    filter(UIN == i) %>%
    pull(CDR3) %>%
    unique()
  ppd <- read.csv(paste0("data/Summary_PPD-reactive-CDR3s_exp-thr3_beta.csv")) %>%
    filter(UIN == i) %>%
    pull(CDR3) %>%
    unique()
  ppd.new <- setdiff(ppd,nil)
  write.table(ppd.new,paste0("data/",i,"_beta_PPD_gr3_cleaned.csv"), sep = ",", row.names = F, col.names = F)
}

### Step 2: down-sampled repertoires ####
### Step 2a: quantify metaclones in blood and TST samples (participants with in vitro data) ####
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(UIN %in% id & tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:4)) {
  print(paste0("Expansion threshold: ",f))
  # make empty result data frame
  result <- data.frame()
  # select TCRs
  dat <- data.b %>%
    filter(duplicate_count >f) %>%
    ungroup()
  # load metaclonotypist reference (mhc II)
  mc <- read.csv("data/TableS4.csv")
  mc <- mc %>% dplyr::select(regex,Vs,index)
  # make empty summary dataframe
  summary <- data.frame()
  # loop through all regex
  for (i in (1:nrow(mc))){
    # define search patterns and metaclone index
    regex <- mc[i,1]
    V <- mc[i,2]
    index <- mc[i,3]
    # look for matches in combined data
    print(paste0("checking metaclone index ",i," of ",nrow(mc)))
    match <- dat %>% filter(str_detect(bioidentity, regex) & str_detect(bioidentity, V))
    # add columns to indicate match and metaclone index
    if (dim(match)[1] != 0){
      match$mc.index <- index
      match$mc.match <- 1
      # add to summary file
      summary <- rbind(summary,match)
    }
  }
  # write to file
  write.csv(summary,paste0("data/metaclonotypist_results_down-sampled_beta_expanded_gr",f,"_n=12.csv"),row.names = F)
}

### Step 2b: metaclones vs. private PPD-CDR3s in Day 2 TST ####
spl <- meta %>%
  filter(UIN %in% id & tissue == "TST_D2") %>%
  pull(sample)

dat <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(sample %in% spl)

# loop through different clone sizes
for (f in c(0:4)){
  summary.ppd <- data.frame()
  # loop through each sample
  for (i in unique(dat$sample)){
    uin <- str_split(i,pattern="_")[[1]][1]
    # select repertoire and calculate total TCR count
    d <- dat %>% filter(sample == i & duplicate_count > f)
    total <- sum(d$duplicate_count)
    # load private ppd reactive sequences
    ppd <- read.csv(paste0("data/",uin,"_beta_PPD_gr3_cleaned.csv"),header=F) %>% pull(V1)
    # quantify %TCRs matching private cdr3s
    ppd.tst <- d %>% filter(junction_aa %in% ppd)
    sum.ppd <- sum(ppd.tst$duplicate_count)
    pct.ppd <- sum.ppd/total*100
    # load metaclone CDR3s
    mc.res <- read.csv(paste0("data/metaclonotypist_results_down-sampled_beta_expanded_gr",f,"_n=12.csv")) %>%
      filter(sample == i)
    # quantify %TCRs matching metaclones
    sum.mc <- sum(mc.res$duplicate_count)
    pct.mc <- sum.mc/total*100
    # only private
    ppd.only <- setdiff(ppd.tst$junction_aa,mc.res$junction_aa)
    dat.ppd.only <- d %>% filter(junction_aa %in% ppd.only)
    sum.ppd.only <- sum(dat.ppd.only$duplicate_count)
    pct.ppd.only <- sum.ppd.only/total*100
    # add to summary
    res <- c(i,uin,total,sum.mc,pct.mc,sum.ppd.only,pct.ppd.only)
    summary.ppd <- rbind(summary.ppd,res)
  }
  colnames(summary.ppd) <- c("Sample","UIN","total.CDR3s","number.mc","pct.mc","number.PPD.only","pct.PPD.only")
  summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")] <- sapply(summary.ppd[c("total.CDR3s","number.mc","pct.mc","number.PPD.only","pct.PPD.only")],as.numeric)
  write.csv(summary.ppd,paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr",f,".csv"),row.names = F)
}

### Step 2c: metaclones vs. private PPD-CDR3s in Day 7 TST ####
spl <- meta %>% filter(UIN %in% id & tissue == "TST_D7") %>% pull(sample)

dat <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(sample %in% spl)

# loop through different clone sizes
for (f in c(0:4)){
  summary.ppd <- data.frame()
  # loop through each sample
  for (i in unique(dat$sample)){
    uin <- str_split(i,pattern="_")[[1]][1]
    # select repertoire and calculate total TCR count
    d <- dat %>% filter(sample == i & duplicate_count > f)
    total <- sum(d$duplicate_count)
    # load private ppd reactive sequences
    ppd <- read.csv(paste0("data/",uin,"_beta_PPD_gr3_cleaned.csv"),header=F) %>% pull(V1)
    # quantify % of all private cdr3s
    ppd.tst <- d %>% filter(junction_aa %in% ppd)
    sum.ppd <- sum(ppd.tst$cdr3_count)
    pct.ppd <- sum.ppd/total*100
    # load metaclone CDR3s
    mc.res <- read.csv(paste0("data/metaclonotypist_results_down-sampled_beta_expanded_gr",f,"_n=12.csv")) %>%
      filter(sample == i)
    # quantify %TCRs matching metaclones
    sum.mc <- sum(mc.res$duplicate_count)
    pct.mc <- sum.mc/total*100
    # only private
    ppd.only <- setdiff(ppd.tst$junction_aa,mc.res$junction_aa)
    dat.ppd.only <- d %>% filter(junction_aa %in% ppd.only)
    sum.ppd.only <- sum(dat.ppd.only$duplicate_count)
    pct.ppd.only <- sum.ppd.only/total*100
    # add to summary
    res <- c(i,uin,total,sum.mc,pct.mc,sum.ppd.only,pct.ppd.only)
    summary.ppd <- rbind(summary.ppd,res)
  }
  colnames(summary.ppd) <- c("Sample","UIN","total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")
  summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")] <- sapply(summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")],as.numeric)
  write.csv(summary.ppd,paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr",f,".csv"),row.names = F)
}

### Step 3: full repertoires ####
### Step 3a: quantify metaclones in blood and TST samples (participants with in vitro data) ####
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(UIN %in% id & tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:4)) {
  print(paste0("Expansion threshold: ",f))
  # make empty result data frame
  result <- data.frame()
  # select TCRs
  dat <- data.b %>%
    filter(duplicate_count >f)
  
  # load metaclonotypist reference
  mc <- read.csv("data/TableS4.csv")
  mc <- mc %>% dplyr::select(regex,Vs,index)
  # make empty summary dataframe
  summary <- data.frame()
  # loop through all regex
  for (i in (1:nrow(mc))){
    # define search patterns and metaclone index
    regex <- mc[i,1]
    V <- mc[i,2]
    index <- mc[i,3]
    # look for matches in combined data
    print(paste0("checking metaclone index ",i," of ",nrow(mc)))
    match <- dat %>% filter(str_detect(bioidentity, regex) & str_detect(bioidentity, V))
    # add columns to indicate match and metaclone index
    if (dim(match)[1] != 0){
      match$mc.index <- index
      match$mc.match <- 1
      # add to summary file
      summary <- rbind(summary,match)
    }
  }
  # write to file
  write.csv(summary,paste0("data/metaclonotypist_results_full-repertoires_beta_expanded_gr",f,"_n=12.csv"),row.names = F)
}

### Step 3b: metaclones vs. private PPD-CDR3s in Day 2 TST ####
spl <- meta %>% filter(UIN %in% id & tissue == "TST_D2") %>% pull(sample)

dat <- fread("data/combined_beta.csv.gz") %>%
  filter(sample %in% spl)

# loop through different clone sizes
for (f in c(0:4)){
  summary.ppd <- data.frame()
  # loop through each sample
  for (i in unique(dat$sample)){
    uin <- str_split(i,pattern="_")[[1]][1]
    # select repertoire and calculate total TCRs
    d <- dat %>% filter(sample == i & duplicate_count > f)
    total <- sum(d$duplicate_count)
    # load private ppd reactive sequences
    ppd <- read.csv(paste0("data/",uin,"_beta_PPD_gr3_cleaned.csv"),header=F) %>% pull(V1)
    # quantify % of all private cdr3s
    ppd.tst <- d %>% filter(junction_aa %in% ppd)
    sum.ppd <- sum(ppd.tst$duplicate_count)
    pct.ppd <- sum.ppd/total*100
    # load metaclone CDR3s
    mc.res <- read.csv(paste0("data/metaclonotypist_results_full-repertoires_beta_expanded_gr",f,"_n=12.csv")) %>%
      filter(sample == i)
    # quantify %TCRs matching metaclones
    sum.mc <- sum(mc.res$duplicate_count)
    pct.mc <- sum.mc/total*100
    # only private
    ppd.only <- setdiff(ppd.tst$junction_aa,mc.res$junction_aa)
    dat.ppd.only <- d %>% filter(junction_aa %in% ppd.only)
    sum.ppd.only <- sum(dat.ppd.only$duplicate_count)
    pct.ppd.only <- sum.ppd.only/total*100
    # add to summary
    res <- c(i,uin,total,sum.mc,pct.mc,sum.ppd.only,pct.ppd.only)
    summary.ppd <- rbind(summary.ppd,res)
  }
  colnames(summary.ppd) <- c("Sample","UIN","total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")
  summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")] <- sapply(summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")],as.numeric)
  write.csv(summary.ppd,paste0("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr",f,".csv"),row.names = F)
}

### Step 3c: metaclones vs. private PPD-CDR3s in Day 7 TST ####
spl <- meta %>% filter(UIN %in% id & tissue == "TST_D7") %>% pull(sample)

dat <- fread("data/combined_beta.csv.gz") %>%
  filter(sample %in% spl)

# loop through different clone sizes
for (f in c(0:4)){
  summary.ppd <- data.frame()
  # loop through each sample
  for (i in unique(dat$sample)){
    uin <- str_split(i,pattern="_")[[1]][1]
    # select repertoire and calculate total TCRs
    d <- dat %>% filter(sample == i & duplicate_count > f)
    total <- sum(d$duplicate_count)
    # load private ppd reactive sequences
    ppd <- read.csv(paste0("data/",uin,"_beta_PPD_gr3_cleaned.csv"),header=F) %>% pull(V1)
    # quantify % of all private cdr3s
    ppd.tst <- d %>% filter(junction_aa %in% ppd)
    sum.ppd <- sum(ppd.tst$duplicate_count)
    pct.ppd <- sum.ppd/total*100
    # load metaclone CDR3s
    mc.res <- read.csv(paste0("data/metaclonotypist_results_full-repertoires_beta_expanded_gr",f,"_n=12.csv")) %>%
      filter(sample == i)
    # quantify %TCRs matching metaclones
    sum.mc <- sum(mc.res$duplicate_count)
    pct.mc <- sum.mc/total*100
    # only private
    ppd.only <- setdiff(ppd.tst$junction_aa,mc.res$junction_aa)
    dat.ppd.only <- d %>% filter(junction_aa %in% ppd.only)
    sum.ppd.only <- sum(dat.ppd.only$duplicate_count)
    pct.ppd.only <- sum.ppd.only/total*100
    # add to summary
    res <- c(i,uin,total,sum.mc,pct.mc,sum.ppd.only,pct.ppd.only)
    summary.ppd <- rbind(summary.ppd,res)
  }
  colnames(summary.ppd) <- c("Sample","UIN","total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")
  summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")] <- sapply(summary.ppd[c("total.TCRs","number.mc","pct.mc","number.PPD.only","pct.PPD.only")],as.numeric)
  write.csv(summary.ppd,paste0("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr",f,".csv"),row.names = F)
}
