# Script_7 (metaclones): tidy Gliph2 output
# - add 'pattern' to each cluster;
# - remove any 'local' patterns;
# - convert pattern to regex;
# - retain only most significant HLA association

library(tidyverse)

# change here
mhc_class <- 'II' # select from 'I' or 'II'

# load significant gliph clusters
gliph <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_output.csv")) 
length(unique(gliph$cluster))

# annotate with 'pattern'
hla <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_HLA.csv"))
clust <- read.csv(paste0("gliph2_beta_mhc",mhc_class,"_cluster.csv"))

clust <- clust %>% 
  select(index,TcRb) %>%
  unique()
hla <- hla %>%
  left_join(clust,relationship = "many-to-many") %>%
  mutate(pattern = ifelse(pattern == "single",TcRb,pattern)) %>%
  select(index,pattern) %>%
  rename(cluster = "index") %>%
  unique()

gliph <- gliph %>%
  left_join(hla, relationship = "many-to-many")

# remove local patterns (if any)
gliph <- gliph %>%
  filter(!str_detect(pattern, "^l"))

# convert pattern into regex
gliph$pattern <- gsub("^g","",gliph$pattern)
gliph$pattern <- gsub("%",".",gliph$pattern)

# keep only most significant association per cluster
gliph$X <- NULL
gliph <- gliph %>% group_by(cluster) %>% slice_min(pvalue) %>% ungroup()

table <- if(mhc_class == "II"){"Table_S6"}else{"Table_S9"}
write.csv(gliph,paste0(table,"_gliph2_beta_mhc",mhc_class,".csv"),row.names = F)
