# Script_7 (metaclones): finalise Gliph2 results
# - convert pattern to regular expression
# - replace 'single' pattern (= fully public clusters) with CDR3

library(tidyverse)

# change here
mhc_class <- 'II' # select from 'I' or 'II'

# load significant gliph clusters
gliph <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_output2.csv")) 

# convert pattern into regex
gliph$pattern <- gsub("%",".",gliph$pattern)

# replace single pattern with CDR3
gliph <- gliph %>%
  mutate(pattern = ifelse(pattern == "single", 
                          paste(unique(unlist(str_split(test1$CDR3s, "\\|")))), 
                          pattern))
  
  
table <- if(mhc_class == "II"){"Table_S6"}else{"Table_S7"}
write.csv(gliph,paste0(table,"_gliph2_beta_mhc",mhc_class,".csv"),row.names = F)