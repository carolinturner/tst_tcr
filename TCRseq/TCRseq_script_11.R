# Abundance of different sets of TCRs in validation datasets

# Validation datasets:
# - PBMC bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1
# - Mtb-reactive T-cells sc-TCRseq: download Supplementary Table 2 from Musvosvi et al. (doi: 10.1038/s41591-022-02110-9)
# - SARS-CoV2-reactive T-cells sc-TCRseq: request access to data from Lindeboom et al. (doi: 10.1038/s41586-024-07575-x)
# - TB Lung sc-TCRseq: download data from GSE253828, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Cancer Lung sc-TCRseq: download data from GSE154826, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Lung/Blood/CD4-T bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1


library(data.table)
library(tidyverse)
library(rstatix)

# Step 0: define different sets of TCRs ####
# all TCRs used for metaclonotypist/Gliph2 analysis
tst.dat <- fread("data/combined_subsampled_5000_10000_beta.csv.gz",na.strings = c("","NA"))
hla.uin <- read.csv("data/TableS3.csv") %>% pull(UIN) %>% unique()
tst.cdr3 <- tst.dat %>% 
  filter(UIN %in% hla.uin & nchar(CDR3B) >5 & clonal_count >1) %>%
  pull(CDR3B) %>%
  unique()

# METACLONOTYPIST
mc <- read.csv("data/TableS4.csv") %>% select(regex,Vs,index)

# GLIPH2 pattern
gliph.pattern <- read.csv("data/TableS6.csv") %>%
  filter(!pattern == "G.T") %>%
  pull(pattern)

# GLIPH2 cdr3
gliph.cdr3 <- read.csv("data/TableS6.csv") %>% pull(CDR3s)

# Step 1: quantify different sets of TCRs in different datasets ####

### DATASET 1: bulk-TCRseq of PPD- and TT-stimulated PBMC ####
dat <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("PBMC_PPD","PBMC_TT"))

## DISCOVERY CDR3s
summary <- dat %>% filter(junction_aa %in% tst.cdr3)
# sum total TCR counts for each tissue
df <- dat %>%
  group_by(tissue) %>% 
  summarise(total.count = sum(duplicate_count))
# sum up total metaclone count for each tissue
mc.dat <- summary %>%
  group_by(tissue) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculate mc.percentages
merge <- inner_join(df,mc.dat) # samples with no metaclones discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_PBMC-bulk_tst-cdr3_results_beta.csv"),row.names = F)

## METACLONOTYPIST
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity, regex) & str_detect(bioidentity, V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total metaclone count for each tissue
mc.dat <- summary %>%
  group_by(tissue) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculate mc.percentages
merge <- inner_join(df,mc.dat) # samples with no metaclones discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_PBMC-bulk_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search pattern
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(bioidentity, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total gliph count for each tissue
gliph.dat <- summary %>%
  group_by(tissue) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculte mc.percentages
merge <- inner_join(df,gliph.dat) # samples with no gliph patterns discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_PBMC-bulk_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 CDR3s
summary <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search pattern
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(bioidentity, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total gliph count for each tissue
gliph.dat <- summary %>%
  group_by(tissue) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculte mc.percentages
merge <- inner_join(df,gliph.dat) # samples with no gliph patterns discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_PBMC-bulk_gliph2-cdr3_results_beta.csv"),row.names = F)

### DATASET 2: sc-TCRseq of Mtb-reactive T cells (Musvosvi, Nat Med, 2023) ####
dat <- read.csv("data/Musvosvi_TableS2.csv",row.names = 1) %>%
  filter(flag == "GOOD") %>%
  select(batch_bc_well,donorId,Cell.Type,Vb,Jb,CDR3b) %>%
  na.omit() %>%
  mutate(bioidentity = paste0(Vb,CDR3b,Jb)) 

## DISCOVERY CDR3s
results <- dat %>% filter(CDR3b %in% tst.cdr3)
results.df <- data.frame(tissue = "Tcells_Mtb",
                         total.count = length(unique(dat$batch_bc_well)),
                         total.mc.count = length(unique(results$batch_bc_well)),
                         mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_tst-cdr3_results_beta.csv"),row.names = F)

## Metaclonotypist
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity,regex) & str_detect(bioidentity,V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# calculate total number of TCRs and metaclones
results <- summary %>%
  select(batch_bc_well,mc.match,mc.index) %>%
  na.omit()

results.df <- data.frame(tissue = "Tcells_Mtb",
                total.count = length(unique(dat$batch_bc_well)),
                total.mc.count = length(unique(results$batch_bc_well)),
                mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search pattern
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(bioidentity, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# calculate total number of TCRs and metaclones
results <- summary %>% select(batch_bc_well,mc.match,mc.index) %>% na.omit()

results.df <- data.frame(tissue = "Tcells_Mtb",
                         total.count = length(unique(dat$batch_bc_well)),
                         total.mc.count = length(unique(results$batch_bc_well)),
                         mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 CDR3s
summary <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search pattern
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(bioidentity, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# calculate total number of TCRs and metaclones
results <- summary %>% select(batch_bc_well,mc.match,mc.index) %>% na.omit()

results.df <- data.frame(tissue = "Tcells_Mtb",
                         total.count = length(unique(dat$batch_bc_well)),
                         total.mc.count = length(unique(results$batch_bc_well)),
                         mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_gliph2-cdr3_results_beta.csv"),row.names = F)

### DATASET 3: sc-TCRseq of SARS-CoV2-reactive T cells (Lindeboom, Nature, 2024) ####
dat <- read.table("data/Lindeboom_df_tcr_pbmcs.txt",header = T, na.strings=c("","NA")) %>%
  filter(has_ir == "True") %>%
  select(barcode,sample_id,patient_id,time_point,covid_status,cell_type,
         IR_VDJ_1_v_gene,IR_VDJ_1_cdr3,IR_VDJ_1_j_gene,IR_VDJ_2_v_gene,IR_VDJ_2_cdr3,IR_VDJ_2_j_gene) %>%
  filter(!if_all(c(IR_VDJ_1_cdr3,IR_VDJ_2_cdr3), is.na)) %>%
  mutate(bioidentity1 = paste0(IR_VDJ_1_v_gene,
                               IR_VDJ_1_cdr3,
                               IR_VDJ_1_j_gene),
         bioidentity2 = paste0(IR_VDJ_2_v_gene, 
                               IR_VDJ_2_cdr3,
                               IR_VDJ_2_j_gene)) # note: up to two beta TCRs per cell

## DISCOVERY CDR3s
results <- dat %>% filter(IR_VDJ_1_cdr3 %in% tst.cdr3 | IR_VDJ_2_cdr3 %in% tst.cdr3)
results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_tst-cdr3_results_beta.csv"),row.names = F)

## METACLONOTYPIST
# make empty summary dataframe 1 (for first beta TCR)
summary1 <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity1, regex) & str_detect(bioidentity1, V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index1 <- index
    match$mc.match1 <- 1
    # add to summary file
    summary1 <- rbind(summary1,match)
  }
}

# make empty summary dataframe 2 (for second beta TCR)
summary2 <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity2, regex) & str_detect(bioidentity2, V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index2 <- index
    match$mc.match2 <- 1
    # add to summary file
    summary2 <- rbind(summary2,match)
  }
}

# merge summaries
results <- full_join(summary1,summary2)
  
# calculate total number of TCRs and metaclones
results <- results %>%
  select(barcode,mc.match1,mc.match2,mc.index1,mc.index2) %>%
  filter(mc.match1 == 1 | mc.match2 ==1)

results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
# make empty summary dataframe 1 (for first TCR)
summary1 <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search patterns and metaclone index
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(IR_VDJ_1_cdr3, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index1 <- motif
    match$mc.match1 <- 1
    # add to summary file
    summary1 <- rbind(summary1,match)
  }
}
# make empty summary dataframe 2 (for second TCR)
summary2 <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search patterns and metaclone index
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(IR_VDJ_2_cdr3, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index2 <- motif
    match$mc.match2 <- 1
    # add to summary file
    summary2 <- rbind(summary2,match)
  }
}
# merge summaries
results <- full_join(summary1,summary2)

# calculate total number of TCRs and metaclones
results <- results %>%
  select(barcode,mc.match1,mc.match2,mc.index1,mc.index2) %>%
  filter(mc.match1 == 1 | mc.match2 ==1)

results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 cdr3s
# make empty summary dataframe 1 (for first TCR)
summary1 <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search patterns and metaclone index
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(IR_VDJ_1_cdr3, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index1 <- motif
    match$mc.match1 <- 1
    # add to summary file
    summary1 <- rbind(summary1,match)
  }
}
# make empty summary dataframe 2 (for second TCR)
summary2 <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search patterns and metaclone index
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(IR_VDJ_2_cdr3, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index2 <- motif
    match$mc.match2 <- 1
    # add to summary file
    summary2 <- rbind(summary2,match)
  }
}
# merge summaries
results <- full_join(summary1,summary2)

# calculate total number of TCRs and metaclones
results <- results %>%
  select(barcode,mc.match1,mc.match2,mc.index1,mc.index2) %>%
  filter(mc.match1 == 1 | mc.match2 ==1)

results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph2-cdr3_results_beta.csv"),row.names = F)

### DATASET 4: sc-TCRseq of TB lung (GSE253828) ####
dat <- read.csv("data/TB_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true" & chain == "TRB") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

## DISCOVERY CDR3s
results <- dat %>% filter(cdr3 %in% tst.cdr3)
results.df <- data.frame(tissue = "scLung_TB",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_TB-lung-sc_tst-cdr3_results_beta.csv"),row.names = F)

## METACLONOTYPIST
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclone index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity,regex) & str_detect(bioidentity,V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()

results.df <- data.frame(tissue = "scLung_TB",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_TB-lung-sc_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search patterns and metaclone index
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(bioidentity,motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}

# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()
results.df <- data.frame(tissue = "scLung_TB",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_TB-lung-sc_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 cdr3
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.cdr3))){
  # define search patterns and metaclone index
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(bioidentity,motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}

# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()
results.df <- data.frame(tissue = "scLung_TB",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_TB-lung-sc_gliph2-cdr3_results_beta.csv"),row.names = F)

### DATASET 5: sc-TCRseq of Cancer lung (GSE154826) ####
dat <- read.csv("data/Cancer_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true" & chain == "TRB") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

## DISCOVERY CDR3s
results <- dat %>% filter(cdr3 %in% tst.cdr3)
results.df <- data.frame(tissue = "scLung_Cancer",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_tst-cdr3_results_beta.csv"),row.names = F)

## METACLONOTYPIST
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclone index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity,regex) & str_detect(bioidentity,V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()

results.df <- data.frame(tissue = "scLung_Cancer",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search patterns and metaclone index
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(bioidentity,motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}

# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()
results.df <- data.frame(tissue = "scLung_Cancer",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 cdr3
summary <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search patterns and metaclone index
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(bioidentity,motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}

# calculate total number of TCRs and metaclones
results <- summary %>% select(barcode,mc.match,mc.index) %>% na.omit()
results.df <- data.frame(tissue = "scLung_Cancer",
                         total.count = length(unique(dat$barcode)),
                         total.mc.count = length(unique(results$barcode)),
                         mc.pct.total = round((length(unique(results$barcode)) / length(unique(dat$barcode)) * 100),3))
# save to file
write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_gliph2-cdr3_results_beta.csv"),row.names = F)

### DATASET 6: bulk-TCRseq of Lung/Blood/CD4-T from TB and cancer ####

# load data and harmonise nomenclature
dat <- fread("data/ImmunoSeq_combined_beta.csv.gz") %>%
  mutate(bioident = gsub("TCRBV","TRBV",bioidentity),
         bioident = gsub("0([1-9])","\\1",bioident))
  
## DISCOVERY CDR3s
summary <- dat %>% filter(aminoAcid %in% tst.cdr3)

# sum total TCR counts for each dataset/group combination
df <- dat %>%
  group_by(dataset, group) %>% 
  summarise(total.count = sum(duplicate_count))
# sum up total metaclone count for each dataset/group combination
mc.dat <- summary %>%
  group_by(dataset, group) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculate mc.percentages
merge <- inner_join(df,mc.dat) # samples with no metaclones discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_tst-cdr3_results_beta.csv"),row.names = F)

## METACLONOTYPIST
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc))){
  # define search patterns and metaclone index
  regex <- mc[i,1]
  V <- mc[i,2]
  index <- mc[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioident, regex) & str_detect(bioident, V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total metaclone count for each dataset/group combination
mc.dat <- summary %>%
  group_by(dataset, group) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculate mc.percentages
merge <- inner_join(df,mc.dat) # samples with no metaclones discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_metaclonotypist_results_beta.csv"),row.names = F)

## GLIPH2 pattern
summary <- data.frame()
# loop through all regex
for (i in (1:length(gliph.pattern))){
  # define search pattern
  motif <- gliph.pattern[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.pattern)))
  match <- dat %>% filter(str_detect(bioident, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total gliph count for each dataset/group combination
gliph.dat <- summary %>%
  group_by(dataset, group) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculte mc.percentages
merge <- inner_join(df,gliph.dat) # samples with no gliph patterns discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_gliph2-pattern_results_beta.csv"),row.names = F)

## GLIPH2 cdr3
# make empty summary dataframe
summary <- data.frame()
# loop through all CDR3 sets
for (i in (1:length(gliph.cdr3))){
  # define search pattern
  motif <- gliph.cdr3[i]
  # look for matches in combined data
  print(paste0("checking gliph index ",i," of ",length(gliph.cdr3)))
  match <- dat %>% filter(str_detect(bioident, motif))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- motif
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# sum up total gliph count for each dataset/group combination
gliph.dat <- summary %>%
  group_by(dataset, group) %>%
  summarise(total.mc.count = sum(duplicate_count))
# merge dataframes and calculte mc.percentages
merge <- inner_join(df,gliph.dat) # samples with no gliph patterns discarded
merge <- merge %>%
  mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
# save to file
write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_gliph2-cdr3_results_beta.csv"),row.names = F)

# Step 2: extract data for plotting and calculate odds ratios ####

# tst-cdr3 results
d1 <- read.csv("data/Validation_PBMC-bulk_tst-cdr3_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_tst-cdr3_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_tst-cdr3_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_tst-cdr3_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_tst-cdr3_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_tst-cdr3_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
tst.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Discovery TST CDR3s") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

# metaclonotypist results
d1 <- read.csv("data/Validation_PBMC-bulk_metaclonotypist_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_metaclonotypist_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_metaclonotypist_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_metaclonotypist_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_metaclonotypist_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_metaclonotypist_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
mc.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Metaclonotypist regex") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 
  
# gliph2-pattern results
d1 <- read.csv("data/Validation_PBMC-bulk_gliph2-pattern_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph2-pattern_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph2-pattern_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph2-pattern_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph2-pattern_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph2-pattern_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph1.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2 pattern") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

# gliph2-pattern results
d1 <- read.csv("data/Validation_PBMC-bulk_gliph2-cdr3_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph2-cdr3_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph2-cdr3_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph2-cdr3_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph2-cdr3_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph2-cdr3_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph2.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2-matching CDR3s") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

## odds ratio calculations
# combine data frames and reformat
df <- rbind(tst.dat,mc.dat,gliph1.dat,gliph2.dat) %>%
  mutate(non.mc.count = total.count - total.mc.count,
         mc.count = total.mc.count) %>%
  select(algorithm,dataset,Stimulant,mc.count,non.mc.count) %>%
  mutate(Stimulant = recode(Stimulant,
                            "PPD" = "TB",
                            "TT" = "nonTB",
                            "Mtb" = "TB",
                            "SARS-CoV2" = "nonTB",
                            "Cancer" = "nonTB",
                            "Lung" = "TB",
                            "PBMC" = "nonTB"))

# odds ratio tst-cdr3
tst.wide <- df %>%
  filter(algorithm == "Discovery TST CDR3s") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(tst.wide,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
tst.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Discovery TST CDR3s") %>%
  rownames_to_column("dataset")

# odds ratio metaclonotypist
mc.wide <- df %>%
  filter(algorithm == "Metaclonotypist regex") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(mc.wide,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
mc.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Metaclonotypist regex") %>%
  rownames_to_column("dataset")

# odds ratio gliph2 pattern
gliph1.wide <- df %>%
  filter(algorithm == "Gliph2 pattern") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph1.wide,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph1.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2 pattern") %>%
  rownames_to_column("dataset")

# odds ratio gliph2 pattern
gliph2.wide <- df %>%
  filter(algorithm == "Gliph2-matching CDR3s") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph2.wide,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph2.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2-matching CDR3s") %>%
  rownames_to_column("dataset")

# save to file
or <- rbind(tst.res,mc.res,gliph1.res,gliph2.res)
write.csv(or,"data/Figure5A_odds-ratio.csv",row.names = F)
