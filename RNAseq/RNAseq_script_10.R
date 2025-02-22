## RNAseq script 9: Module analysis

# download Table S1 from paper and save in working directory as 'modules.csv'

library(tidyverse)

tpm <- read.csv("tpm_PC0.001_log2_genesymbol_dedup.csv",row.names = 1)
modules <- read.csv("modules.csv") 

# calculate module score per sample
summary <- data.frame()

for (i in 1:nrow(modules)){
  mod <- modules[i,1]
  genes <- modules[i,2] %>% str_split(",") %>% unlist()
  score <- tpm %>%
    filter(row.names(tpm) %in% genes) %>%
    summarise(across(all_of(names(tpm)),mean))
  summary <- rbind(summary,score)
}
row.names(summary) <- modules$module

summary.t <- as.data.frame(t(summary)) %>% rownames_to_column("sample")
write.csv(summary.t,"Modules_summary.csv",row.names = F)

# scale module scores using saline samples as control
meta <- read.csv("RNAseq_metadata.csv") %>% select(-UIN)
dat <- left_join(summary.t,meta)  
  
# convert into long format
d <- dat %>%
  pivot_longer(cols = c(-Stimulant,-sample),names_to = "module_name", values_to = "moduleTPM") %>% 
  mutate_at(c("moduleTPM"), as.numeric)

# filter control values (=Saline) and calculate mean and standard deviation
c <- filter(d, Stimulant == "saline") %>% 
  group_by(module_name) %>%
  summarise(mean(moduleTPM), sd(moduleTPM))
colnames(c) <- c("module_name", "mean.moduleTPM", "sd.moduleTPM")

# merge with original df by module_name and calculate Z-scores
d <- left_join(d, c, by = "module_name") %>%
  mutate(module_Zscore = (moduleTPM-mean.moduleTPM)/sd.moduleTPM)
write.csv(d, "Modules_Z-scores.csv",row.names = F)
