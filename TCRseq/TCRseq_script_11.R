# Metaclone and gliph pattern abundance in validation datasets
# - PBMC bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1
# - Mtb-reactive T-cells sc-TCRseq: download Supplementary Table 2 from Musvosvi et al. (doi: 10.1038/s41591-022-02110-9)
# - SARS-CoV2-reactive T-cells sc-TCRseq: request access to data from Lindeboom et al. (doi: 10.1038/s41586-024-07575-x)
# - TB Lung sc-TCRseq: download data from GSE253828, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Cancer Lung sc-TCRseq: download data from GSE154826, process with CellRanger's vdj pipeline and combine 'filtered_contig_annotations.csv' files from all samples into one file
# - Lung/Blood/CD4-T bulk-TCRseq: download from UCL RDR and process with TCRseq_script_1

library(data.table)
library(tidyverse)
library(rstatix)

# Step 1: quantify metaclones in different datasets using metaclonotypist or gliph2 motifs ####

# dataset 1: bulk-TCRseq of PPD- and TT-stimulated PBMC ####
dat.a <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("PBMC_PPD","PBMC_TT"))
dat.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("PBMC_PPD","PBMC_TT"))

# search for metaclones
for (chain in c("alpha","beta")){
  print(paste0("Chain: ",chain))
  # assign dataset
  dat <- if(chain == "alpha"){dat.a}else{dat.b}
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
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
  
  # sum total TCR counts for each tissue
  dat <- dat %>%
    group_by(tissue) %>% 
    summarise(total.count = sum(duplicate_count))
  # sum up total metaclone count for each tissue
  mc.dat <- summary %>%
    group_by(tissue) %>%
    summarise(total.mc.count = sum(duplicate_count))
  # merge dataframes and calculate mc.percentages
  merge <- inner_join(dat,mc.dat) # samples with no metaclones discarded
  merge <- merge %>%
    mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
  # save to file
  write.csv(merge,paste0("data/Validation_PBMC-bulk_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search pattern
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat.b %>% filter(str_detect(bioidentity, regex.beta))
      # add columns to indicate match and metaclone index
      if (dim(match)[1] != 0){
        match$mc.index <- regex[i]
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
    merge <- inner_join(dat,gliph.dat) # samples with no gliph patterns discarded
    merge <- merge %>%
      mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
    # save to file
    write.csv(merge,paste0("data/Validation_PBMC-bulk_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# dataset 2: sc-TCRseq of Mtb-reactive T cells (Musvosvi, Nat Med, 2023) ####
dat <- read.csv("data/Musvosvi_TableS2.csv",row.names = 1) %>%
  filter(flag == "GOOD") %>%
  select(batch_bc_well,donorId,Cell.Type,Vb,Jb,CDR3b,Va,Ja,CDR3a,alt.Va,alt.Ja,alt.CDR3a)

# split by chain and add column with merged V_gene and CDR3 annotation
dat.alpha <- dat %>%
  select(-Vb,-Jb,-CDR3b) %>%
  mutate(bioidentity = paste0(Va,CDR3a,Ja),
         bioidentity.alt = paste0(alt.Va,alt.CDR3a,alt.Ja)) # note: up to two alpha TCRs per cell
dat.beta <- dat %>%
  select(-Va,-Ja,-CDR3a,-alt.Va,-alt.Ja,-alt.CDR3a) %>%
  mutate(bioidentity = paste0(Vb,CDR3b,Jb))

# search for metaclones
for (chain in c("alpha","beta")){
  print(paste0("Chain: ",chain))
  # assign dataset
  dat <- if(chain == "alpha"){dat.alpha}else{dat.beta}
  
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
  mc <- mc %>% dplyr::select(regex,Vs,index)
  # make empty summary dataframe 1
  summary1 <- data.frame()
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
      summary1 <- rbind(summary1,match)
    }
  }
  if(chain == "alpha"){
    # make empty summary dataframe 2
    summary2 <- data.frame()
    # loop through all regex
    for (i in (1:nrow(mc))){
      # define search patterns and metaclone index
      regex <- mc[i,1]
      V <- mc[i,2]
      index <- mc[i,3]
      # look for matches in combined data
      print(paste0("checking metaclone index ",i," of ",nrow(mc)))
      match <- dat %>% filter(str_detect(bioidentity.alt, regex) & str_detect(bioidentity.alt, V))
      # add columns to indicate match and metaclone index
      if (dim(match)[1] != 0){
        match$mc.index.alt <- index
        match$mc.match.alt <- 1
        # add to summary file
        summary2 <- rbind(summary2,match)
      }
    }
  }
  # merge summaries
  results <- if(chain == "alpha"){full_join(summary1,summary2)}else{summary1}
  
  # calculate total number of TCRs and metaclones
  results <- if(chain == "alpha"){
    results %>%
      select(batch_bc_well,mc.match,mc.index,mc.match.alt,mc.index.alt) %>%
      filter(!if_all(c(mc.match,mc.match.alt),is.na))
  } else {
    results %>% select(batch_bc_well,mc.match,mc.index) %>%
      na.omit()
  }
  
  d <- if(chain == "alpha"){
    dat %>% filter(!if_all(c(CDR3a,alt.CDR3a),is.na))
  } else {
    dat %>% na.omit()
  }
  
  results.df <- data.frame(tissue = "Tcells_Mtb",
                  total.count = length(unique(d$batch_bc_well)),
                  total.mc.count = length(unique(results$batch_bc_well)),
                  mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(d$batch_bc_well)) * 100),3))
  # save to file
  write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search pattern
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat %>% filter(str_detect(bioidentity, regex.beta))
      # add columns to indicate match and metaclone index
      if (dim(match)[1] != 0){
        match$mc.index <- regex[i]
        match$mc.match <- 1
        # add to summary file
        summary <- rbind(summary,match)
      }
    }
    # calculate total number of TCRs and metaclones
    results <- summary %>% select(batch_bc_well,mc.match,mc.index) %>% na.omit()
    d <- dat %>% na.omit()
    
    results.df <- data.frame(tissue = "Tcells_Mtb",
                             total.count = length(unique(d$batch_bc_well)),
                             total.mc.count = length(unique(results$batch_bc_well)),
                             mc.pct.total = round((length(unique(results$batch_bc_well)) / length(unique(d$batch_bc_well)) * 100),3))
    # save to file
    write.csv(results.df,paste0("data/Validation_Mtb-Tcells-sc_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# dataset 3: sc-TCRseq of SARS-CoV2-reactive T cells (Lindeboom, Nature, 2024) ####
dat <- read.table("data/Lindeboom_df_tcr_pbmcs.txt",header = T, na.strings=c("","NA")) %>%
  filter(has_ir == "True") %>%
  select(barcode,sample_id,patient_id,time_point,covid_status,cell_type,
         IR_VDJ_1_v_gene,IR_VDJ_1_cdr3,IR_VDJ_1_j_gene,IR_VDJ_2_v_gene,IR_VDJ_2_cdr3,IR_VDJ_2_j_gene,
         IR_VJ_1_v_gene,IR_VJ_1_cdr3,IR_VJ_1_j_gene,IR_VJ_2_v_gene,IR_VJ_2_cdr3,IR_VJ_2_j_gene)

# split by chain and add column with merged V_gene and CDR3 annotation
dat.alpha <- dat %>%
  select(-IR_VDJ_1_cdr3,-IR_VDJ_1_v_gene,-IR_VDJ_1_j_gene,-IR_VDJ_2_cdr3,-IR_VDJ_2_j_gene,-IR_VDJ_2_v_gene) %>%
  mutate(bioidentity1 = paste0(IR_VJ_1_v_gene,
                               IR_VJ_1_cdr3,
                               IR_VJ_1_j_gene),
         bioidentity2 = paste0(IR_VJ_2_v_gene, 
                               IR_VJ_2_cdr3,
                               IR_VJ_2_j_gene)) # note: up to two alpha TCRs per cell
dat.beta <- dat %>%
  select(-IR_VJ_1_v_gene,-IR_VJ_1_cdr3,-IR_VJ_1_j_gene,-IR_VJ_2_cdr3,-IR_VJ_2_j_gene,-IR_VJ_2_v_gene) %>%
  mutate(bioidentity1 = paste0(IR_VDJ_1_v_gene,
                               IR_VDJ_1_cdr3,
                               IR_VDJ_1_j_gene),
         bioidentity2 = paste0(IR_VDJ_2_v_gene, 
                               IR_VDJ_2_cdr3,
                               IR_VDJ_2_j_gene)) # note: up to two beta TCRs per cell

# search for metaclones
for (chain in c("alpha","beta")){
  print(paste0("Chain: ",chain))
  # assign dataset
  dat <- if(chain == "alpha"){dat.alpha}else{dat.beta}
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
  mc <- mc %>% dplyr::select(regex,Vs,index)
  # make empty summary dataframe 1 (for first alpha or beta TCR)
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
  
  # make empty summary dataframe 2 (for second alpha or beta TCR)
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
  
  d <- if(chain == "alpha"){
    dat %>% filter(!if_all(c(IR_VJ_1_cdr3,IR_VJ_2_cdr3),is.na))
  } else {
    dat %>% filter(!if_all(c(IR_VDJ_1_cdr3,IR_VDJ_2_cdr3), is.na))
  }
  
  results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                           total.count = length(unique(d$barcode)),
                           total.mc.count = length(unique(results$barcode)),
                           mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
  # save to file
  write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe 1 (for first TCR)
    summary1 <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search patterns and metaclone index
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat %>% filter(str_detect(IR_VDJ_1_cdr3, regex.beta))
      # add columns to indicate match and metaclone index
      if (dim(match)[1] != 0){
        match$mc.index1 <- index
        match$mc.match1 <- 1
        # add to summary file
        summary1 <- rbind(summary1,match)
      }
    }
    # make empty summary dataframe 2 (for second TCR)
    summary2 <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search patterns and metaclone index
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat %>% filter(str_detect(IR_VDJ_2_cdr3, regex.beta))
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
    
    d <- dat %>% filter(!if_all(c(IR_VDJ_1_cdr3,IR_VDJ_2_cdr3), is.na))
    
    results.df <- data.frame(tissue = "Tcells_SARS-CoV2",
                             total.count = length(unique(d$barcode)),
                             total.mc.count = length(unique(results$barcode)),
                             mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
    # save to file
    write.csv(results.df,paste0("data/Validation_SARS-CoV2-Tcells-sc_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# dataset 4: sc-TCRseq of TB lung (GSE253828) ####
dat <- read.csv("data/TB_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

for (chain in c("alpha","beta")){
  print(paste0("Testing chain: ",chain))
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
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
  
  d <- if(chain == "alpha"){
    dat %>% filter(chain == "TRA")
  } else {
    dat %>% filter(chain == "TRB")
  }
  
  results.df <- data.frame(tissue = "scLung_TB",
                           total.count = length(unique(d$barcode)),
                           total.mc.count = length(unique(results$barcode)),
                           mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
  # save to file
  write.csv(results.df,paste0("data/Validation_TB-lung-sc_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search patterns and metaclone index
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat %>% filter(str_detect(bioidentity,regex.beta))
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
    d <- dat %>% filter(chain == "TRB")
    results.df <- data.frame(tissue = "scLung_TB",
                             total.count = length(unique(d$barcode)),
                             total.mc.count = length(unique(results$barcode)),
                             mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
    # save to file
    write.csv(results.df,paste0("data/Validation_TB-lung-sc_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# dataset 5: sc-TCRseq of Cancer lung (GSE154826) ####
dat <- read.csv("data/Cancer_lung_scTCRseq_combined.csv") %>%
  filter(is_cell == "true") %>%
  mutate(bioidentity = paste0(v_gene,cdr3,j_gene))

for (chain in c("alpha","beta")){
  print(paste0("Testing chain: ",chain))
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
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
  
  d <- if(chain == "alpha"){
    dat %>% filter(chain == "TRA")
  } else {
    dat %>% filter(chain == "TRB")
  }
  
  results.df <- data.frame(tissue = "scLung_Cancer",
                           total.count = length(unique(d$barcode)),
                           total.mc.count = length(unique(results$barcode)),
                           mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
  # save to file
  write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search patterns and metaclone index
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat %>% filter(str_detect(bioidentity,regex.beta))
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
    d <- dat %>% filter(chain == "TRB")
    results.df <- data.frame(tissue = "scLung_Cancer",
                             total.count = length(unique(d$barcode)),
                             total.mc.count = length(unique(results$barcode)),
                             mc.pct.total = round((length(unique(results$barcode)) / length(unique(d$barcode)) * 100),3))
    # save to file
    write.csv(results.df,paste0("data/Validation_Cancer-lung-sc_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# dataset 6: bulk-TCRseq of Lung/Blood/CD4-T from TB and cancer ####

# load data and harmonise nomenclature
dat.a <- fread("data/ImmunoSeq_combined_alpha.csv.gz") %>%
  mutate(bioident = gsub("TCRAV","TRAV",bioidentity),
         bioident = gsub("0([1-9])","\\1",bioident))
dat.b <- fread("data/ImmunoSeq_combined_beta.csv.gz") %>%
  mutate(bioident = gsub("TCRBV","TRBV",bioidentity),
         bioident = gsub("0([1-9])","\\1",bioident))
         
# search for metaclones
for (chain in c("alpha","beta")){
  print(paste0("Chain: ",chain))
  # assign dataset
  dat <- if(chain == "alpha"){dat.a}else{dat.b}
  # load metaclonotypist reference
  mc <- if(chain == "alpha"){
    read.csv("data/TableS5.csv")
  } else {
    read.csv("data/TableS4.csv")
  }
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
  
  # sum total TCR counts for each dataset/group combination
  dat <- dat %>%
    group_by(dataset, group) %>% 
    summarise(total.count = sum(duplicate_count))
  # sum up total metaclone count for each dataset/group combination
  mc.dat <- summary %>%
    group_by(dataset, group) %>%
    summarise(total.mc.count = sum(duplicate_count))
  # merge dataframes and calculate mc.percentages
  merge <- inner_join(dat,mc.dat) # samples with no metaclones discarded
  merge <- merge %>%
    mutate(mc.pct.total = round(total.mc.count/total.count*100,3))
  # save to file
  write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_metaclonotypist_results_",chain,".csv"),row.names = F)
  
  # gliph2
  if(chain == "beta"){
    # load gliph2 reference and remove 'G.T' pattern (as too unspecific)
    gliph <- read.csv("data/TableS6.csv")
    regex <- gliph$pattern
    regex <- regex[!regex == "G.T"]
    # make empty summary dataframe
    summary <- data.frame()
    # loop through all regex
    for (i in (1:length(regex))){
      # define search pattern
      regex.beta <- regex[i]
      # look for matches in combined data
      print(paste0("checking gliph index ",i," of ",length(regex)))
      match <- dat.b %>% filter(str_detect(bioident, regex.beta))
      # add columns to indicate match and metaclone index
      if (dim(match)[1] != 0){
        match$mc.index <- regex[i]
        match$mc.match <- 1
        # add to summary file
        summary <- rbind(summary,match)
      }
    }
    
    # sum up total gliph count for each dataset/group combination
    gliph.dat <- summary %>%
      group_by(dataset, group) %>%
      summarise(total.gliph.count = sum(duplicate_count))
    # merge dataframes and calculte mc.percentages
    merge <- inner_join(dat,gliph.dat) # samples with no gliph patterns discarded
    merge <- merge %>%
      mutate(gliph.mc.total = round(total.mc.count/total.count*100,3))
    # save to file
    write.csv(merge,paste0("data/Validation_ImmunoSeq-bulk_gliph2_results_",chain,".csv"),row.names = F)
  }
}

# Step 2: extract data for plotting and calculate odds ratios ####

# beta ####
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
  mutate(algorithm = "Metaclonotypist") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 
  
# gliph2 results
d1 <- read.csv("data/Validation_PBMC-bulk_gliph2_results_beta.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_gliph2_results_beta.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_gliph2_results_beta.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_gliph2_results_beta.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_gliph2_results_beta.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_gliph2_results_beta.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
gliph.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Gliph2") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) 

## data for plotting
dat <- rbind(mc.dat,gliph.dat) %>%
  select(dataset,Stimulant,mc.pct.total,algorithm) %>%
  mutate(dataset = recode(dataset,
                          'PBMC' = "PBMC (bulk-TCRseq)",
                          'Tcells' = "T-cells (sc-TCRseq)",
                          'scLung' = "Lung (sc-TCRseq)",
                          'Lung' = "Lung (bulk-TCRseq)",
                          'Blood' = "Blood (bulk-TCRseq)",
                          'CD4-T' = "CD4-T (bulk-TCRseq)"),
         Stimulant = recode(Stimulant,
                            'Lung' = "TB lung",
                            'PBMC' = "TB blood"))
write.csv(dat,"data/Figure5A_percentage.csv",row.names = F)

## odds ratio calculations
# combine data frames and reformat
df <- rbind(mc.dat,gliph.dat) %>%
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

# odds ratio metaclonotypist
mc.wide <- df %>%
  filter(algorithm == "Metaclonotypist") %>%
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
  mutate(algorithm = "Metaclonotypist") %>%
  rownames_to_column("dataset")

# odds ratio gliph2
gliph.wide <- df %>%
  filter(algorithm == "Gliph2") %>%
  pivot_wider(names_from = Stimulant, values_from = c(mc.count,non.mc.count)) %>%
  as.data.frame() %>%
  column_to_rownames("dataset") %>%
  select(mc.count_TB,non.mc.count_TB,mc.count_nonTB,non.mc.count_nonTB)

# do fisher test
res <- t(apply(gliph.wide,
               1,
               function(x){
                 x1 <- fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T))
                 c(x1$estimate,x1$conf.int,x1$p.value)}
))
colnames(res) <- c("oddsratio","CI_lower","CI_higher","pval")
gliph.res <- as.data.frame(res) %>%
  adjust_pvalue(method = "fdr") %>% # correct for multiple testing
  mutate_at(vars(oddsratio, CI_lower, CI_higher), ~round(., 1)) %>%
  mutate(algorithm = "Gliph2") %>%
  rownames_to_column("dataset")

# save to file
or <- rbind(mc.res,gliph.res)
write.csv(or,"data/Figure5A_odds-ratio.csv",row.names = F)

# alpha ####
d1 <- read.csv("data/Validation_PBMC-bulk_metaclonotypist_results_alpha.csv")
d2 <- read.csv("data/Validation_Mtb-Tcells-sc_metaclonotypist_results_alpha.csv")
d3 <- read.csv("data/Validation_SARS-CoV2-Tcells-sc_metaclonotypist_results_alpha.csv")
d4 <- read.csv("data/Validation_TB-lung-sc_metaclonotypist_results_alpha.csv")
d5 <- read.csv("data/Validation_Cancer-lung-sc_metaclonotypist_results_alpha.csv")
d6 <- read.csv("data/Validation_ImmunoSeq-bulk_metaclonotypist_results_alpha.csv") %>%
  mutate(tissue = paste0(dataset,"_",group)) %>%
  select(-dataset,-group)
mc.dat <- rbind(d1,d2,d3,d4,d5,d6) %>%
  mutate(algorithm = "Metaclonotypist") %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant"))

## data for plotting
dat <- mc.dat %>%
  select(tissue,mc.pct.total) %>%
  separate_wider_delim(tissue, "_", names = c("dataset","Stimulant")) %>%
  mutate(dataset = recode(dataset,
                          'PBMC' = "PBMC (bulk-TCRseq)",
                          'Tcells' = "T-cells (sc-TCRseq)",
                          'scLung' = "Lung (sc-TCRseq)",
                          'Lung' = "Lung (bulk-TCRseq)",
                          'Blood' = "Blood (bulk-TCRseq)"))
write.csv(dat,"data/FigureS8A_percentage.csv",row.names = F)
## odds ratio calculations
# reformat data frame
df <- mc.dat %>%
  mutate(non.mc.count = total.count - total.mc.count,
         mc.count = total.mc.count) %>%
  select(algorithm,dataset,Stimulant,mc.count,non.mc.count) %>%
  mutate(Stimulant = recode(Stimulant,
                            "PPD" = "TB",
                            "TT" = "nonTB",
                            "Mtb" = "TB",
                            "SARS-CoV2" = "nonTB",
                            "Cancer" = "nonTB"))

# odds ratio metaclonotypist
mc.wide <- df %>%
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
  mutate(algorithm = "Metaclonotypist") %>%
  rownames_to_column("dataset")

# save to file
write.csv(mc.res,"data/FigureS8A_odds-ratio.csv",row.names = F)