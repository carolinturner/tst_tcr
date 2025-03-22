# Script_7 (metaclones): add regex pattern to Gliph2 clusters

library(tidyverse)

# change here
mhc_class <- 'II' # select from 'I' or 'II'

# load significant gliph clusters
gliph <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_output2.csv")) 
length(unique(gliph$cluster))

# annotate with 'pattern'
hla <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_HLA.csv")) %>%
  rename(cluster = index,
         pattern_type = pattern) %>%
  select(cluster,pattern_type) %>%
  unique()

gliph <- gliph %>%
  left_join(hla) %>%
  select(-pattern_type,-X) %>%
  arrange(cluster) %>%
  rownames_to_column("index")

# convert pattern into regex
gliph$pattern <- gsub("%",".",gliph$pattern)

table <- if(mhc_class == "II"){"Table_S6"}else{"Table_S7"}
write.csv(gliph,paste0(table,"_gliph2_beta_mhc",mhc_class,".csv"),row.names = F)