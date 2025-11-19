## TCRseq_script_10: Effect of filtering on sensitivity and specificity of Gliph2 clusters

library(tidyverse)
library(data.table)
library(rstatix)

### Step 1: define sets of CDR3s ####
d1 <- read.csv("data/gliph2_filter1.csv") # Gliph2 + Filter 1
d2 <- read.csv("data/gliph2_filter1+2.csv") %>% filter(pvalue < 0.05) # Gliph2 + Filter 1 + Filter 2 (p<0.05) 
d3 <- read.csv("data/gliph2_filter1+2.csv") %>% filter(significant == "True") # Gliph2 + Filter 1 + Filter 2 (FDR<0.1)
d4 <- read.csv("data/FileS7.csv") # Gliph2 + Filter 2 (FDR<0.1)
# in contrast to Gliph2 output shown in Figure 5A, pattern 'G.T' is retained in the analysis here

# summary stats to manually populate table for Figure S10A
# number of clusters
d1.clust <- length(unique(d1$index)) # 1837
d2.clust <- length(unique(d2$cluster)) # 1342
d3.clust <- length(unique(d3$cluster)) # 106
d4.clust <- length(unique(d4$index)) # 128

# number of CDR3s
d1.cdr3 <- unique(d1$TcRb) # 6546
d2.cdr3 <- separate_longer_delim(d2,CDR3s,"|") %>% pull(CDR3s) %>% unique() # 5022
d3.cdr3 <- separate_longer_delim(d3,CDR3s,"|") %>% pull(CDR3s) %>% unique() # 476
d4.cdr3 <- separate_longer_delim(d4,CDR3s,"|") %>% pull(CDR3s) %>% unique() # 476

# overlap between 476 CDR3s retained with Filter 2 +/- Filter 1 
all(d3.cdr3 %in% d4.cdr3) # FALSE
length(intersect(d3.cdr3,d4.cdr3)) # 401

### Step 2: find matches in independent validation datasets ####

# Validation datasets:
# - PBMC bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1
# - Mtb-reactive T-cells sc-TCRseq: download Supplementary Table 2 from Musvosvi et al. (doi: 10.1038/s41591-022-02110-9)
# - SARS-CoV2-reactive T-cells sc-TCRseq: request access to data from Lindeboom et al. (doi: 10.1038/s41586-024-07575-x)
# - TB Lung sc-TCRseq: download data from GSE253828, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Cancer Lung sc-TCRseq: download data from GSE154826, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Lung/Blood/CD4-T bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1

## DATASET 1: bulk-TCRseq of PPD- and TT-stimulated PBMC ####
dat <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("PBMC_PPD","PBMC_TT"))

# find matching cdr3s
summary1 <- dat %>% filter(junction_aa %in% d1.cdr3)
summary2 <- dat %>% filter(junction_aa %in% d2.cdr3)
summary3 <- dat %>% filter(junction_aa %in% d3.cdr3)
summary4 <- dat %>% filter(junction_aa %in% d4.cdr3)

# sum total TCR counts for each tissue
df <- dat %>% group_by(tissue) %>% summarise(total.count = sum(duplicate_count))
# sum up total metaclone count for each tissue
mc.dat1 <- summary1 %>% group_by(tissue) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat2 <- summary2 %>% group_by(tissue) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat3 <- summary3 %>% group_by(tissue) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat4 <- summary4 %>% group_by(tissue) %>% summarise(total.mc.count = sum(duplicate_count))

# merge dataframes, calculate mc.percentages and save to file
merge1 <- inner_join(df,mc.dat1) # samples with no metaclones discarded
merge1 <- merge1 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge1,paste0("data/Validation_PBMC-bulk_gliph-filter_set1_results_beta.csv"),row.names = F)

merge2 <- inner_join(df,mc.dat2) # samples with no metaclones discarded
merge2 <- merge2 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge2,paste0("data/Validation_PBMC-bulk_gliph-filter_set2_results_beta.csv"),row.names = F)

merge3 <- inner_join(df,mc.dat3) # samples with no metaclones discarded
merge3 <- merge3 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge3,paste0("data/Validation_PBMC-bulk_gliph-filter_set3_results_beta.csv"),row.names = F)

merge4 <- inner_join(df,mc.dat4) # samples with no metaclones discarded
merge4 <- merge4 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge4,paste0("data/Validation_PBMC-bulk_gliph-filter_set4_results_beta.csv"),row.names = F)

## DATASET 2: sc-TCRseq of Mtb-reactive T cells (Musvosvi, Nat Med, 2023) ####
dat <- read.csv("data/Musvosvi_TableS2.csv",row.names = 1) %>%
  filter(flag == "GOOD") %>%
  select(batch_bc_well,donorId,Cell.Type,Vb,Jb,CDR3b) %>%
  na.omit() %>%
  mutate(bioidentity = paste0(Vb,CDR3b,Jb)) 

# find matching cdr3s
summary1 <- dat %>% filter(CDR3b %in% d1.cdr3)
summary2 <- dat %>% filter(CDR3b %in% d2.cdr3)
summary3 <- dat %>% filter(CDR3b %in% d3.cdr3)
summary4 <- dat %>% filter(CDR3b %in% d4.cdr3)

# calculate mc.percentage and save to file
results1 <- data.frame(tissue = "Tcells_Mtb",
                       total.count = length(unique(dat$batch_bc_well)),
                       total.mc.count = length(unique(summary1$batch_bc_well)),
                       mc.pct.total = round((length(unique(summary1$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
write.csv(results1,paste0("data/Validation_Mtb-Tcells-sc_gliph-filter_set1_results_beta.csv"),row.names = F)

results2 <- data.frame(tissue = "Tcells_Mtb",
                       total.count = length(unique(dat$batch_bc_well)),
                       total.mc.count = length(unique(summary2$batch_bc_well)),
                       mc.pct.total = round((length(unique(summary2$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
write.csv(results2,paste0("data/Validation_Mtb-Tcells-sc_gliph-filter_set2_results_beta.csv"),row.names = F)

results3 <- data.frame(tissue = "Tcells_Mtb",
                       total.count = length(unique(dat$batch_bc_well)),
                       total.mc.count = length(unique(summary3$batch_bc_well)),
                       mc.pct.total = round((length(unique(summary3$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
write.csv(results3,paste0("data/Validation_Mtb-Tcells-sc_gliph-filter_set3_results_beta.csv"),row.names = F)

results4 <- data.frame(tissue = "Tcells_Mtb",
                       total.count = length(unique(dat$batch_bc_well)),
                       total.mc.count = length(unique(summary4$batch_bc_well)),
                       mc.pct.total = round((length(unique(summary4$batch_bc_well)) / length(unique(dat$batch_bc_well)) * 100),3))
write.csv(results4,paste0("data/Validation_Mtb-Tcells-sc_gliph-filter_set4_results_beta.csv"),row.names = F)

## DATASET 3: sc-TCRseq of SARS-CoV2-reactive T cells (Lindeboom, Nature, 2024) ####
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

# find matching cdr3s
summary1 <- dat %>% filter(IR_VDJ_1_cdr3 %in% d1.cdr3 | IR_VDJ_2_cdr3 %in% d1.cdr3)
summary2 <- dat %>% filter(IR_VDJ_1_cdr3 %in% d2.cdr3 | IR_VDJ_2_cdr3 %in% d2.cdr3)
summary3 <- dat %>% filter(IR_VDJ_1_cdr3 %in% d3.cdr3 | IR_VDJ_2_cdr3 %in% d3.cdr3)
summary4 <- dat %>% filter(IR_VDJ_1_cdr3 %in% d4.cdr3 | IR_VDJ_2_cdr3 %in% d4.cdr3)

# calculate mc.percentage and save to file
results1 <- data.frame(tissue = "Tcells_SARS-CoV2",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary1$barcode)),
                       mc.pct.total = round((length(unique(summary1$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results1,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set1_results_beta.csv"),row.names = F)

results2 <- data.frame(tissue = "Tcells_SARS-CoV2",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary2$barcode)),
                       mc.pct.total = round((length(unique(summary2$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results2,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set2_results_beta.csv"),row.names = F)

results3 <- data.frame(tissue = "Tcells_SARS-CoV2",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary3$barcode)),
                       mc.pct.total = round((length(unique(summary3$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results3,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set3_results_beta.csv"),row.names = F)

results4 <- data.frame(tissue = "Tcells_SARS-CoV2",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary4$barcode)),
                       mc.pct.total = round((length(unique(summary4$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results4,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set4_results_beta.csv"),row.names = F)

## DATASET 4: sc-TCRseq of TB lung (GSE253828) ####
dat <- read.csv("data/TB_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true" & chain == "TRB") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

# find matching cdr3s
summary1 <- dat %>% filter(cdr3 %in% d1.cdr3)
summary2 <- dat %>% filter(cdr3 %in% d2.cdr3)
summary3 <- dat %>% filter(cdr3 %in% d3.cdr3)
summary4 <- dat %>% filter(cdr3 %in% d4.cdr3)

# calculate mc.percentage and save to file
results1 <- data.frame(tissue = "scLung_TB",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary1$barcode)),
                       mc.pct.total = round((length(unique(summary1$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results1,paste0("data/Validation_TB-lung-sc_gliph-filter_set1_results_beta.csv"),row.names = F)

results2 <- data.frame(tissue = "scLung_TB",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary2$barcode)),
                       mc.pct.total = round((length(unique(summary2$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results2,paste0("data/Validation_TB-lung-sc_gliph-filter_set2_results_beta.csv"),row.names = F)

results3 <- data.frame(tissue = "scLung_TB",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary3$barcode)),
                       mc.pct.total = round((length(unique(summary3$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results3,paste0("data/Validation_TB-lung-sc_gliph-filter_set3_results_beta.csv"),row.names = F)

results4 <- data.frame(tissue = "scLung_TB",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary4$barcode)),
                       mc.pct.total = round((length(unique(summary4$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results4,paste0("data/Validation_TB-lung-sc_gliph-filter_set4_results_beta.csv"),row.names = F)

## DATASET 5: sc-TCRseq of Cancer lung (GSE154826) ####
dat <- read.csv("data/Cancer_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true" & chain == "TRB") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

# find matching cdr3s
summary1 <- dat %>% filter(cdr3 %in% d1.cdr3)
summary2 <- dat %>% filter(cdr3 %in% d2.cdr3)
summary3 <- dat %>% filter(cdr3 %in% d3.cdr3)
summary4 <- dat %>% filter(cdr3 %in% d4.cdr3)

# calculate mc.percentage and save to file
results1 <- data.frame(tissue = "scLung_Cancer",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary1$barcode)),
                       mc.pct.total = round((length(unique(summary1$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results1,paste0("data/Validation_Cancer-lung-sc_gliph-filter_set1_results_beta.csv"),row.names = F)

results2 <- data.frame(tissue = "scLung_Cancer",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary2$barcode)),
                       mc.pct.total = round((length(unique(summary2$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results2,paste0("data/Validation_Cancer-lung-sc_gliph-filter_set2_results_beta.csv"),row.names = F)

results3 <- data.frame(tissue = "scLung_Cancer",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary3$barcode)),
                       mc.pct.total = round((length(unique(summary3$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results3,paste0("data/Validation_Cancer-lung-sc_gliph-filter_set3_results_beta.csv"),row.names = F)

results4 <- data.frame(tissue = "scLung_Cancer",
                       total.count = length(unique(dat$barcode)),
                       total.mc.count = length(unique(summary4$barcode)),
                       mc.pct.total = round((length(unique(summary4$barcode)) / length(unique(dat$barcode)) * 100),3))
write.csv(results4,paste0("data/Validation_Cancer-lung-sc_gliph-filter_set4_results_beta.csv"),row.names = F)

## DATASET 6: bulk-TCRseq of Lung/Blood/CD4-T from TB and cancer ####

# load data and harmonise nomenclature
dat <- fread("data/ImmunoSeq_combined_beta.csv.gz") %>%
  mutate(bioident = gsub("TCRBV","TRBV",bioidentity),
         bioident = gsub("0([1-9])","\\1",bioident))

# find matching cdr3s
summary1 <- dat %>% filter(aminoAcid %in% d1.cdr3)
summary2 <- dat %>% filter(aminoAcid %in% d2.cdr3)
summary3 <- dat %>% filter(aminoAcid %in% d3.cdr3)
summary4 <- dat %>% filter(aminoAcid %in% d4.cdr3)

# sum total TCR counts for each dataset/group combination
df <- dat %>% group_by(dataset, group) %>% summarise(total.count = sum(duplicate_count))
# sum up total metaclone count for each dataset/group combination
mc.dat1 <- summary1 %>% group_by(dataset, group) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat2 <- summary2 %>% group_by(dataset, group) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat3 <- summary3 %>% group_by(dataset, group) %>% summarise(total.mc.count = sum(duplicate_count))
mc.dat4 <- summary4 %>% group_by(dataset, group) %>% summarise(total.mc.count = sum(duplicate_count))

# merge dataframes, calculate mc.percentages and save to file
merge1 <- inner_join(df,mc.dat1) # samples with no metaclones discarded
merge1 <- merge1 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge1,paste0("data/Validation_ImmunoSeq-bulk_gliph-filter_set1_results_beta.csv"),row.names = F)

merge2 <- inner_join(df,mc.dat2) # samples with no metaclones discarded
merge2 <- merge2 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge2,paste0("data/Validation_ImmunoSeq-bulk_gliph-filter_set2_results_beta.csv"),row.names = F)

merge3 <- inner_join(df,mc.dat3) # samples with no metaclones discarded
merge3 <- merge3 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge3,paste0("data/Validation_ImmunoSeq-bulk_gliph-filter_set3_results_beta.csv"),row.names = F)

merge4 <- inner_join(df,mc.dat4) # samples with no metaclones discarded
merge4 <- merge4 %>% mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
write.csv(merge4,paste0("data/Validation_ImmunoSeq-bulk_gliph-filter_set4_results_beta.csv"),row.names = F)

### Step 3: extract data for plotting and calculate odds ratios ####

# results set 1
d1 <- read.csv("data/Validation_PBMC-bulk_gliph-filter_set1_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph-filter_set1_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set1_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph-filter_set1_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph-filter_set1_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph-filter_set1_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph1.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2 + Filter 1") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

# results set 2
d1 <- read.csv("data/Validation_PBMC-bulk_gliph-filter_set2_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph-filter_set2_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set2_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph-filter_set2_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph-filter_set2_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph-filter_set2_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph2.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2 + Filter 1 + Filter 2 (p<0.05)") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

# results set 3
d1 <- read.csv("data/Validation_PBMC-bulk_gliph-filter_set3_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph-filter_set3_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set3_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph-filter_set3_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph-filter_set3_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph-filter_set3_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph3.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2 + Filter 1 + Filter 2 (FDR<0.1)") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

# results set 4
d1 <- read.csv("data/Validation_PBMC-bulk_gliph-filter_set4_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph-filter_set4_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph-filter_set4_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph-filter_set4_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph-filter_set4_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph-filter_set4_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph4.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2 + Filter 2 (FDR<0.1)") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 


## odds ratio calculations
# combine data frames and reformat
df <- rbind(gliph1.dat,gliph2.dat,gliph3.dat,gliph4.dat) %>%
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

# odds ratio set1
gliph1 <- df %>%
  filter(algorithm == "Gliph2 + Filter 1") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph1,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph1.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2 + Filter 1") %>%
  rownames_to_column("dataset")

# odds ratio set2
gliph2 <- df %>%
  filter(algorithm == "Gliph2 + Filter 1 + Filter 2 (p<0.05)") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph2,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph2.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2 + Filter 1 + Filter 2 (p<0.05)") %>%
  rownames_to_column("dataset")

# odds ratio set3
gliph3 <- df %>%
  filter(algorithm == "Gliph2 + Filter 1 + Filter 2 (FDR<0.1)") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph3,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph3.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2 + Filter 1 + Filter 2 (FDR<0.1)") %>%
  rownames_to_column("dataset")

# odds ratio set4
gliph4 <- df %>%
  filter(algorithm == "Gliph2 + Filter 2 (FDR<0.1)") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph4,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph4.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2 + Filter 2 (FDR<0.1)") %>%
  rownames_to_column("dataset")

# save to file
or <- rbind(gliph1.res,gliph2.res,gliph3.res,gliph4.res)
write.csv(or,"data/Figure10B_odds-ratio_gliph-filters.csv",row.names = F)
