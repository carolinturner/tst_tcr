## TCRseq_script_16: Assess publicity of Mtb-reactive TCRs in day 7 TST

library(data.table)
library(tidyverse)

# Step 1: full repertoires ####
# load TST_D7 beta chain repertoires 
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue == "TST_D7") 

# Calculate publicity (using different expansion thresholds)
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs & count number of participants
  dat <- data.b %>% filter(duplicate_count >f)
  n <- length(unique(dat$sample))
  
  # load references
  mc.dat <- read.csv(paste0("data/metaclonotypist_results_full-repertoires_beta_expanded_gr",f,".csv")) %>%
    filter(tissue == "TST_D7")
  mtb.dat <- read.csv(paste0("data/Published-Ag-search_full-repertoires_expanded_gr",f,"_beta.csv")) %>%
    filter(Mtb == "Mtb" & tissue == "TST_D7")
  
  # calculate publicity
  cdr3 <- dat %>%
    select(sample,junction_aa) %>%
    unique() %>%
    group_by(junction_aa) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "CDR3") %>%
    rename(Id = "junction_aa") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  mc <- mc.dat %>%
    select(sample,mc.index) %>%
    unique() %>%
    group_by(mc.index) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "metaclone") %>%
    rename(Id = "mc.index") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  mtb <- mtb.dat %>%
    select(sample,junction_aa) %>%
    unique() %>%
    group_by(junction_aa) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "CDR3 with published Mtb reactivity") %>%
    rename(Id = "junction_aa") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  # calculate cumulative proportion of people captured by most common metaclones/cdr3s
  for (i in c(1:nrow(mc))){
    df <- subset(mc.dat, mc.index %in% c(mc[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    mc[i,5] <- prop  
  }
  colnames(mc)[5] <- "cum.prop.people"
  
  for (i in c(1:nrow(mc))){ # only need to match number of metaclones
    df <- subset(dat, junction_aa %in% c(cdr3[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    cdr3[i,5] <- prop  
  }
  colnames(cdr3)[5] <- "cum.prop.people"
  
  for (i in c(1:nrow(mc))){ # only need to match number of metaclones
    df <- subset(mtb.dat, junction_aa %in% c(mtb[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    mtb[i,5] <- prop  
  }
  colnames(mtb)[5] <- "cum.prop.people"
  
  # merge & save
  summary <- rbind(mc,cdr3,mtb) %>% filter(rank <= nrow(mc))
  write.csv(summary,paste0("data/Publicity_all-vs-published-vs-metaclone_full-repertoires_beta_expanded_gr",f,".csv"))
}


# Step 2: down-sampled repertoires ####
# load TST_D7 beta chain repertoires 
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue == "TST_D7")

# Calculate publicity (using different expansion thresholds)
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs & count number of participants
  dat <- data.b %>% filter(duplicate_count >f)
  n <- length(unique(dat$sample))
  
  # load references
  mc.dat <- read.csv(paste0("data/metaclonotypist_results_down-sampled_beta_expanded_gr",f,".csv")) %>%
    filter(tissue == "TST_D7")
  mtb.dat <- read.csv(paste0("data/Published-Ag-search_down-sampled_expanded_gr",f,"_beta.csv")) %>%
    filter(Mtb == "Mtb" & tissue == "TST_D7")
  
  # calculate publicity
  cdr3 <- dat %>%
    select(sample,junction_aa) %>%
    unique() %>%
    group_by(junction_aa) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "CDR3") %>%
    rename(Id = "junction_aa") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  mc <- mc.dat %>%
    select(sample,mc.index) %>%
    unique() %>%
    group_by(mc.index) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "metaclone") %>%
    rename(Id = "mc.index") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  mtb <- mtb.dat %>%
    select(sample,junction_aa) %>%
    unique() %>%
    group_by(junction_aa) %>%
    mutate(publicity = length(unique(sample)),
           dataset = "CDR3 with published Mtb reactivity") %>%
    rename(Id = "junction_aa") %>%
    ungroup() %>%
    select(-sample) %>%
    unique() %>%
    arrange(-publicity) %>%
    rowid_to_column("rank")
  
  # calculate cumulative proportion of people captured by most common metaclones/cdr3s
  for (i in c(1:nrow(mc))){
    df <- subset(mc.dat, mc.index %in% c(mc[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    mc[i,5] <- prop  
  }
  colnames(mc)[5] <- "cum.prop.people"
  
  for (i in c(1:nrow(mc))){ # only need to match number of metaclones
    df <- subset(dat, junction_aa %in% c(cdr3[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    cdr3[i,5] <- prop  
  }
  colnames(cdr3)[5] <- "cum.prop.people"
  
  for (i in c(1:nrow(mc))){ # only need to match number of metaclones
    df <- subset(mtb.dat, junction_aa %in% c(mtb[1:i,2]$Id))
    prop <- (length(unique(df$sample)))/n
    mtb[i,5] <- prop  
  }
  colnames(mtb)[5] <- "cum.prop.people"
  
  # merge & save
  summary <- rbind(mc,cdr3,mtb) %>% filter(rank <= nrow(mc))
  write.csv(summary,paste0("data/Publicity_all-vs-published-vs-metaclone_down-sampled_beta_expanded_gr",f,".csv"))
}
