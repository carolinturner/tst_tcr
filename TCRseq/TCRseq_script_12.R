library(tidyverse)
library(data.table)

# Step 1: extract metaclone CDR3s ####
# beta mhcII
b <- read.csv("data/TableS4.csv")

mc.b <- as.vector(as.matrix(b[,"CDR3s"]))
mc.b <- unlist(strsplit(mc.b,"\\|"))
mc.b <- as.data.frame(unique(mc.b))
colnames(mc.b) <- "CDR3"

write.table(mc.b,"data/beta_metaclone_CDR3s.csv",sep=",",col.names = F,row.names = F,quote = F)

# alpha mhcII
a <- read.csv("data/TableS5.csv")

mc.a <- as.vector(as.matrix(a[,"CDR3s"]))
mc.a <- unlist(strsplit(mc.a,"\\|"))
mc.a <- as.data.frame(unique(mc.a))
colnames(mc.a) <- "CDR3"

write.table(mc.a,"data/alpha_metaclone_CDR3s.csv",sep=",",col.names = F,row.names = F,quote = F)

# Step 2: down-sample published Mtb-reactive CDR3s to same number of metaclone CDR3s ####
mtb.a <- read.csv("data/TableS2.csv") %>%
  filter(reactivity == "Mtb" & chain == "alpha") %>%
  pull(CDR3) %>%
  unique()
mtb.b <- read.csv("data/TableS2.csv") %>%
  filter(reactivity == "Mtb" & chain == "beta") %>%
  pull(CDR3) %>%
  unique()

set.seed(12345)
mtb.sample.a <- sample(mtb.a,nrow(mc.a))
mtb.sample.b <- sample(mtb.b,nrow(mc.b))

# Step 3: quantify abundance in full repertoires ####
data.a <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()
data.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  for (chain in c("alpha","beta")){
    dat <- if(chain == "alpha"){data.a}else{data.b}
    mtb.sample <- if(chain =="alpha"){mtb.sample.a}else{mtb.sample.b}
    mc <- if(chain == "alpha"){mc.a %>% pull(CDR3)}else{mc.b %>% pull(CDR3)}
    # select TCRs
    dat <- dat %>% filter(cdr3_count >f)
    dat$mtb <- ifelse(dat$junction_aa %in% mtb.sample, "Mtb", NA)
    dat$mc <- ifelse(dat$junction_aa %in% mc, "mc", NA)

    # sum total/unique counts for each sample
    d <- dat %>%
      group_by(tissue,sample) %>% 
      summarise(total.count = sum(duplicate_count),
                unique.count = length(duplicate_count)) %>%
      ungroup()
    
    # sum total/unique counts of mc and mtb hits for each sample
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
    write.csv(summary,paste0("data/Summary_percentage_",chain,"_size-matched-mc-mtb_search_full-repertoires_expanded_gr",f,".csv"))
  }
}

# Step 4: quantify abundance in down-sampled repertoires ####
data.a <- fread("data/combined_subsampled_alpha.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()
data.b <- fread("data/combined_subsampled_beta.csv.gz") %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) %>%
  group_by(sample,junction_aa) %>%
  mutate(cdr3_count = sum(duplicate_count)) %>%
  ungroup()

# loop through different expansion thresholds
for (f in c(0:1)) {
  print(paste0("Expansion threshold: ",f))
  for (chain in c("alpha","beta")){
    dat <- if(chain == "alpha"){data.a}else{data.b}
    mtb.sample <- if(chain =="alpha"){mtb.sample.a}else{mtb.sample.b}
    mc <- if(chain == "alpha"){mc.a %>% pull(CDR3)}else{mc.b %>% pull(CDR3)}
    # select TCRs
    dat <- dat %>% filter(cdr3_count >f)
    dat$mtb <- ifelse(dat$junction_aa %in% mtb.sample, "Mtb", NA)
    dat$mc <- ifelse(dat$junction_aa %in% mc, "mc", NA)

    # sum total/unique counts for each sample
    d <- dat %>%
      group_by(tissue,sample) %>% 
      summarise(total.count = sum(duplicate_count),
                unique.count = length(duplicate_count)) %>%
      ungroup()
    
    # sum total/unique counts of mc and mtb hits for each sample
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
    write.csv(summary,paste0("data/Summary_percentage_",chain,"_size-matched-mc-mtb_search_down-sampled_expanded_gr",f,".csv"))
  }
}
