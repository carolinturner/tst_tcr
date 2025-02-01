# TCRseq_script_9: quantify metaclones in blood and TST samples

library(data.table)
library(tidyverse)

# full repertoires ####
data.a <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # make empty result data frame
  result <- data.frame()
  # select TCRs
  dat.a <- data.a %>%
    group_by(sample,junction_aa) %>%
    mutate(count_by_cdr3 = sum(duplicate_count)) %>%
    filter(count_by_cdr3 >f) %>%
    ungroup()
  dat.b <- data.b %>%
    group_by(sample,junction_aa) %>%
    mutate(count_by_cdr3 = sum(duplicate_count)) %>%
    filter(count_by_cdr3 >f) %>%
    ungroup()
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    dat <- if(chain == "alpha"){dat.a}else{dat.b}
    # load metaclonotypist reference
    chain.ref <- if(chain == "alpha"){
      read.csv("data/TableS5.csv")
    }else{
      read.csv("data/TableS4.csv")
    }
    chain.ref <- chain.ref %>% dplyr::select(regex,Vs,index)
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:nrow(chain.ref))){
      # define search patterns and metaclone index
      regex <- chain.ref[i,1]
      V <- chain.ref[i,2]
      index <- chain.ref[i,3]
      # look for matches in combined data
      print(paste0("checking metaclone index ",i," of ",nrow(chain.ref)))
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
    write.csv(summary,paste0("data/metaclonotypist_results_full-repertoires_",chain,"_expanded_gr",f,".csv"),row.names = F)
    
    # calculate percentage of metaclones amongst total CDR3s
    dat <- dat %>%
      group_by(tissue,sample) %>%
      summarise(total.count = sum(duplicate_count)) %>%
      ungroup()
    summary <- summary %>%
      group_by(tissue,sample) %>%
      summarise(mc.count = sum(duplicate_count)) %>%
      ungroup()
    merge <- left_join(dat,summary) %>%
      mutate(mc.pct = mc.count/total.count*100,
             mc.pct = replace_na(mc.pct,0),
             Chain = chain) %>%
      select(Chain,tissue,sample,mc.pct)
    result <- rbind(result,merge)
  }
  write.csv(result,paste0("data/summary_metaclone-abundance_full-repertoires_expanded_gr",f,".csv"),row.names = F)
}

# down-sampled repertoires ####
data.a <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # make empty result data frame
  result <- data.frame()
  # select TCRs
  dat.a <- data.a %>%
    group_by(sample,junction_aa) %>%
    mutate(count_by_cdr3 = sum(duplicate_count)) %>%
    filter(count_by_cdr3 >f) %>%
    ungroup()
  dat.b <- data.b %>%
    group_by(sample,junction_aa) %>%
    mutate(count_by_cdr3 = sum(duplicate_count)) %>%
    filter(count_by_cdr3 >f) %>%
    ungroup()
  
  for (chain in c("alpha","beta")){
    print(paste0("Chain: ",chain))
    # assign dataset
    dat <- if(chain == "alpha"){dat.a}else{dat.b}
    # load metaclonotypist reference
    chain.ref <- if(chain == "alpha"){
      read.csv("data/TableS5.csv")
    }else{
      read.csv("data/TableS4.csv")
    }
    chain.ref <- chain.ref %>% dplyr::select(regex,Vs,index)
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:nrow(chain.ref))){
      # define search patterns and metaclone index
      regex <- chain.ref[i,1]
      V <- chain.ref[i,2]
      index <- chain.ref[i,3]
      # look for matches in combined data
      print(paste0("checking metaclone index ",i," of ",nrow(chain.ref)))
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
    write.csv(summary,paste0("data/metaclonotypist_results_down-sampled_",chain,"_expanded_gr",f,".csv"),row.names = F)
    
    # calculate percentage of metaclones amongst total CDR3s
    dat <- dat %>%
      group_by(tissue,sample) %>%
      summarise(total.count = sum(duplicate_count)) %>%
      ungroup()
    summary <- summary %>%
      group_by(tissue,sample) %>%
      summarise(mc.count = sum(duplicate_count)) %>%
      ungroup()
    merge <- left_join(dat,summary) %>%
      mutate(mc.pct = mc.count/total.count*100,
             mc.pct = replace_na(mc.pct,0),
             Chain = chain) %>%
      select(Chain,tissue,sample,mc.pct)
    result <- rbind(result,merge)
  }
  write.csv(result,paste0("data/summary_metaclone-abundance_down-sampled_expanded_gr",f,".csv"),row.names = F)
}
