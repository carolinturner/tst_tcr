library(tidyverse)
library(data.table)

# Step 1: extract metaclone CDR3s (mhc II) ####
mc2 <- read.csv("data/TableS4.csv")

mc.b <- as.vector(as.matrix(mc2[,"CDR3s"]))
mc.b <- unlist(strsplit(mc.b,"\\|"))
mc.b <- as.data.frame(unique(mc.b))
colnames(mc.b) <- "CDR3"
write.table(mc.b,"data/classII_beta_metaclone_CDR3s.csv",sep=",",col.names = F,row.names = F,quote = F)
mc.cdr3 <- mc.b %>% pull(CDR3)

# Step 2: down-sample published Mtb-reactive CDR3s to same number of metaclone CDR3s ####
mtb.b <- read.csv("data/TableS2.csv") %>%
  filter(reactivity == "Mtb" & chain == "beta") %>%
  pull(CDR3) %>%
  unique()

set.seed(12345)
mtb.sample <- sample(mtb.b,nrow(mc.b))

# Step 3: quantify abundance in full repertoires ####
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) 

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs
  dat <- data.b %>% filter(duplicate_count >f)
  dat$mtb <- ifelse(dat$junction_aa %in% mtb.sample, "Mtb", NA)
  dat$mc <- ifelse(dat$junction_aa %in% mc.cdr3, "mc", NA)
  # sum total/unique counts for each sample
  d <- dat %>%
    group_by(tissue,sample) %>% 
    summarise(total.count = sum(duplicate_count)) %>%
    ungroup()
  
  # sum total counts of mc and mtb hits for each sample
  mtb <- dat %>%
    filter(mtb == "Mtb") %>%
    group_by(tissue,sample) %>%
    summarise(total.mtb.count = sum(duplicate_count)) %>%
    ungroup()
  mc <- dat %>%
    filter(mc == "mc") %>%
    group_by(tissue,sample) %>%
    summarise(total.mc.count = sum(duplicate_count)) %>%
    ungroup()
  
  # merge dataframes and calculate %
  merge.mtb <- inner_join(d,mtb) %>%
    mutate(pct = total.mtb.count/total.count*100,
           TCR = "published") %>%
    select(tissue,sample,TCR,pct)
  merge.mc <- inner_join(d,mc) %>%
    mutate(pct = total.mc.count/total.count*100,
           TCR = "metaclone") %>%
    select(tissue,sample,TCR,pct)
  
  # merge into one summary dataframe and prep for plotting
  summary <- rbind(merge.mtb,merge.mc)
  write.csv(summary,paste0("data/Summary_percentage_beta_size-matched-mc-mtb_search_full-repertoires_expanded_gr",f,".csv"))
}


# Step 4: quantify abundance in down-sampled repertoires ####
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7"))

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs
  dat <- data.b %>% filter(duplicate_count >f)
  dat$mtb <- ifelse(dat$junction_aa %in% mtb.sample, "Mtb", NA)
  dat$mc <- ifelse(dat$junction_aa %in% mc.cdr3, "mc", NA)
    # sum total counts for each sample
  d <- dat %>%
    group_by(tissue,sample) %>% 
    summarise(total.count = sum(duplicate_count)) %>%
    ungroup()
  
  # sum total counts of mc and mtb hits for each sample
  mtb <- dat %>%
    filter(mtb == "Mtb") %>%
    group_by(tissue,sample) %>%
    summarise(total.mtb.count = sum(duplicate_count)) %>%
    ungroup()
  mc <- dat %>%
    filter(mc == "mc") %>%
    group_by(tissue,sample) %>%
    summarise(total.mc.count = sum(duplicate_count)) %>%
    ungroup()
  
  # merge dataframes and calculate %
  merge.mtb <- inner_join(d,mtb) %>%
    mutate(pct = total.mtb.count/total.count*100,
           TCR = "published") %>%
    select(tissue,sample,TCR,pct)
  merge.mc <- inner_join(d,mc) %>%
    mutate(pct = total.mc.count/total.count*100,
           TCR = "metaclone") %>%
    select(tissue,sample,TCR,pct)
  
  # merge into one summary dataframe and prep for plotting
  summary <- rbind(merge.mtb,merge.mc)
  write.csv(summary,paste0("data/Summary_percentage_beta_size-matched-mc-mtb_search_down-sampled_expanded_gr",f,".csv"))
}
