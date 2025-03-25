# TCRseq_script_6: Define private antigen-reactive TCRs

## AIM: Define TCR CDR3 clones that are reactive to in-vitro stimulation with antigen.
# Reactive clones are defined by expansion compared to paired blood sample, in at least one in-vitro replicate.
# Log2 expansion threshold is user-defined (eg. exp=3 selects clones that are expanded 8-fold or more).
# PPD- and TT-reactive CDR3s are cleaned by removing CDR3s that are also reactive to stimulation with Nil in any participant.
# CDR3s that react to both PPD and TT are retained (cross-reactivity).

library(data.table)
library(tidyverse)

# load metadata and identify UINs with in-vitro stimulation data
meta <- read.csv("data/TCRseq_metadata.csv")
id <- meta %>% filter(tissue == "PBMC_PPD") %>% pull(UIN) %>% unique()

# load data and select relevant samples
alpha <- fread("data/combined_alpha.csv.gz") %>%
  filter(UIN %in% id & !tissue %in% c("TST_D2","TST_D7"))
beta <- fread("data/combined_beta.csv.gz") %>%
  filter(UIN %in% id & !tissue %in% c("TST_D2","TST_D7"))

# set expansion threshold
exp <- 3

# Step 1: Define reactive in-vitro CDR3s ####
for(chain in c("alpha","beta")){
  print(paste0("Processing chain: ",chain))
  dat <- if(chain == "alpha"){alpha}else{beta}
  # make empty summary dataframes
  nil.df <- data.frame()
  ppd.df <- data.frame()
  tt.df <- data.frame()
  
  # loop through subjects
  for(i in id){
    print(paste0("Processing UIN: ",i))
    # blood reference data
    blood <- dat %>%
      filter(UIN == i & tissue == "Blood") %>%
      group_by(UIN,junction_aa) %>%
      summarise(total = sum(duplicate_count)) %>%
      ungroup()
    colnames(blood) <- c("UIN","CDR3","Blood")
    # calculate median CDR3 count in blood (usually =1)
    median_blood <- median(blood$Blood)
    # copy dataframe in preparation of merging with PBMC data
    df <- blood
    
    # PBMC data
    pbmc <- dat %>%
      filter(UIN == i & !tissue == "Blood") %>%
      group_by(UIN,sample,replicate,junction_aa) %>%
      summarise(total = sum(duplicate_count)) %>%
      ungroup()
    colnames(pbmc) <- c("UIN","sample","replicate","CDR3","PBMC")

    # merge PBMC data with blood data
    df <- merge(df,pbmc,by=c("CDR3","UIN"), all=TRUE)
    
    # set absent CDR3s in blood to median CDR3 count in blood, then remove CDR3s absent from PBMC data
    df$Blood[is.na(df$Blood)] <- median_blood
    df <- df %>% drop_na()
    
    # log2 transform counts
    df.log <- df %>% mutate(across(where(is.numeric), ~log2(.)))

    # subtract blood counts from PBMC data and filter for CDR3s with stimulation index >= expansion threshold
    df.log <- df.log %>%
      mutate(SI = PBMC - Blood) %>%
      filter(SI >= exp)
    
    # split reactive CDR3s by antigen
    nil <- df.log %>% filter(grepl("Nil",sample))
    ppd <- df.log %>% filter(grepl("PPD",sample))
    tt <- df.log %>% filter(grepl("TT",sample))
    
    # add to summary files
    nil.df <- rbind(nil.df,nil)
    ppd.df <- rbind(ppd.df,ppd)
    tt.df <- rbind(tt.df,tt)
  }
  # save summary files
  write.csv(nil.df, paste0("data/Summary_Nil-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"),row.names = F)
  write.csv(ppd.df, paste0("data/Summary_PPD-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"),row.names = F)
  write.csv(tt.df, paste0("data/Summary_TT-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"),row.names = F)
}

# Step 2: Tidy up and clean lists of reactive CDR3s ####
for (chain in c("alpha","beta")){
  # load summary files
  nil <- read.csv(paste0("data/Summary_Nil-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"))
  ppd <- read.csv(paste0("data/Summary_PPD-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"))
  tt <- read.csv(paste0("data/Summary_TT-reactive-CDR3s_exp-thr",exp,"_",chain,".csv"))
  # retain only unique CDR3s
  nil <- nil %>% pull(CDR3) %>% unique()
  ppd <- ppd %>% pull(CDR3) %>% unique()
  tt <- tt %>% pull(CDR3) %>% unique()
  # clean ppd- and tt-reactive CDR3s by removing CDR3s that also expand in unstimulated cultures
  tt.clean <- setdiff(tt,nil)
  ppd.clean <- setdiff(ppd,nil)
  # save to file
  write.table(ppd.clean,paste0("data/PPD_",chain,".csv"),sep=",",row.names=F,col.names=F)
  write.table(tt.clean,paste0("data/TT_",chain,".csv"),sep=",",row.names=F,col.names=F)
}

# Step 3: Publicness of PPD-reactive CDR3s ####
summary <- data.frame(publicness=factor())

for (chain in c("alpha","beta")){
  # load list of cleaned PPD-reactive CDR3s
  list <- as.character(read.csv(paste0("data/PPD_",chain,".csv"),header=F)$V1)
  # load integrated data and retain only cleaned CDR3s
  dat <- read.csv(paste0("data/Summary_PPD-reactive-CDR3s_exp-thr",exp,"_",chain,".csv")) %>%
    select(UIN,CDR3) %>%
    filter(CDR3 %in% list) 
  # pivot wider
  d <- dat %>%
    group_by(CDR3,UIN) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(expanded = 1) %>%
    pivot_wider(id_cols = CDR3, names_from = UIN, values_from = expanded)
  # add publicness column
  d$publicness <- rowSums(!is.na(d[-1]))
  # calculate frequency distribution of public CDR3s
  freq <- data.frame(table(d$publicness))
  colnames(freq) <- c("publicness",chain)
  # add to summary
  summary <- merge(summary,freq,all=TRUE)
  # make heatmap data by removing publicness column and converting NA to 0
  d$publicness <- NULL
  d[is.na(d)] <- 0
  write.csv(d,paste0("data/Heatmap_data_",chain,".csv"),row.names = F)
}
write.csv(summary, "data/Publicness_in-vitro_PPD.csv")
