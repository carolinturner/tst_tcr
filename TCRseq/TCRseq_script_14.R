library(data.table)
library(tidyverse)

# Step 1: full repertoires ####
# load TST_D7 repertoires and add cdr3_count
data.a <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue == "TST_D7") %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue == "TST_D7") %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()

# Calculate publicity (using different expansion thresholds)
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs & count number of participants
  dat.a <- data.a %>% filter(cdr3_count >f) 
  n.a <- length(unique(dat.a$sample))
  
  dat.b <- data.b %>% filter(cdr3_count >f)
  n.b <- length(unique(dat.b$sample))
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    if(chain == "alpha"){
      dat <- dat.a
      n <- n.a
    }
    if(chain == "beta"){
      dat <- dat.b
      n <- n.b
    }
    # load references
    mc.dat <- read.csv(paste0("data/metaclonotypist_results_full-repertoires_",chain,"_expanded_gr",f,".csv")) %>%
      filter(tissue == "TST_D7")
    mtb.dat <- read.csv(paste0("data/Published-Ag-search_full-repertoires_expanded_gr",f,"_",chain,".csv")) %>%
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
    write.csv(summary,paste0("data/Publicity_all-vs-published-vs-metaclone_full-repertoires_",chain,"_expanded_gr",f,".csv"))
  }
}


# Step 2: down-sampled repertoires ####
# load TST_D7 repertoires and add cdr3_count
data.a <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue == "TST_D7") %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue == "TST_D7") %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()

# Calculate publicity (using different expansion thresholds)
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs & count number of participants
  dat.a <- data.a %>% filter(cdr3_count >f) 
  n.a <- length(unique(dat.a$sample))
  
  dat.b <- data.b %>% filter(cdr3_count >f)
  n.b <- length(unique(dat.b$sample))
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    if(chain == "alpha"){
      dat <- dat.a
      n <- n.a
    }
    if(chain == "beta"){
      dat <- dat.b
      n <- n.b
    }
    # load references
    mc.dat <- read.csv(paste0("data/metaclonotypist_results_down-sampled_",chain,"_expanded_gr",f,".csv")) %>%
      filter(tissue == "TST_D7")
    mtb.dat <- read.csv(paste0("data/Published-Ag-search_down-sampled_expanded_gr",f,"_",chain,".csv")) %>%
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
    write.csv(summary,paste0("data/Publicity_all-vs-published-vs-metaclone_down-sampled_",chain,"_expanded_gr",f,".csv"))
  }
}

