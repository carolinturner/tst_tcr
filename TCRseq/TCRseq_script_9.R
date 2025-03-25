# TCRseq_script_9: quantify class II metaclones in blood and TST samples

library(data.table)
library(tidyverse)

# full repertoire ####
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:1)) {
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
  write.csv(summary,paste0("data/metaclonotypist_results_full-repertoires_beta_expanded_gr",f,".csv"),row.names = F)
  
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
           mc.pct = replace_na(mc.pct,0)) %>%
    select(tissue,sample,mc.pct)
  result <- rbind(result,merge)
  
  write.csv(result,paste0("data/summary_metaclone-abundance_full-repertoires_beta_expanded_gr",f,".csv"),row.names = F)
}

# down-sampled repertoires ####
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:1)) {
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
  write.csv(summary,paste0("data/metaclonotypist_results_down-sampled_beta_expanded_gr",f,".csv"),row.names = F)
  
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
           mc.pct = replace_na(mc.pct,0)) %>%
    select(tissue,sample,mc.pct)
  result <- rbind(result,merge)
  
  write.csv(result,paste0("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr",f,".csv"),row.names = F)
}
