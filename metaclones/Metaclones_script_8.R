## Metaclones_script_8: Check metaclone CDR3s against VDJdb Mtb TCRs
### use files generated here as input to make Venn diagram with deepVenn for Figure 4F
### use output in line 29 to populate the table in Supplementary Figure 9A

library(tidyverse)

# download VDJdb reference at https://github.com/antigenomics/vdjdb-db/releases
## used here: vdjdb-2024-06-13.zip

# load VDJdb reference
ref <- read.delim("data/vdjdb_full.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "") %>%
  filter(mhc.class == "MHCII" & antigen.species %in% c("M.tuberculosis","Mtb")) %>%
  select(cdr3.beta,v.beta,antigen.species,antigen.epitope,antigen.gene)
write.csv(ref,"data/VDJdb_Mtb_MHCII_beta.csv",row.names = F) 

# load gliph2 CDR3s 
mc <- read.csv("data/FileS7.csv")
mc_sep <- separate_longer_delim(mc,CDR3s,"|") %>% unique()
d <- mc_sep %>% select(CDR3s) %>% unique() # 476 unique CDR3s
write.csv(d,"data/gliph2_unique_cdr3s.csv",row.names = F) 

# load metaclonotypist CDR3s 
mc <- read.csv("data/FileS5.csv")
mc_sep <- separate_longer_delim(mc,CDR3s,"|") %>% unique()
d <- mc_sep %>% select(CDR3s) %>% unique() # 1392 unique CDR3s
write.csv(d,"data/metaclonotypist_unique_cdr3s.csv",row.names = F) 

# find matches to manually populate the table in Figure S9A
match <- mc_sep %>% dplyr::rename(cdr3.beta = "CDR3s") %>% inner_join(ref)