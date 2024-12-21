## correlation matrix sc modules vs CCND1_proliferation

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# load module scores
scores <- read.csv("Modules_Summary.csv")

# load metadata
meta <- read.csv("RNAseq_metadata.csv")
meta.stim <- meta %>% select(sample,Stimulant)


### Z-score scaled module scores (using saline samples as control)
# add stimulant column to module score dataframe
scores <- scores %>%
  left_join(meta.stim) %>%
  select(sample,Stimulant,everything())

# pivot longer
d <- scores %>%
  pivot_longer(3:ncol(scores),names_to = "module_name", values_to = "moduleTPM")

# filter control values (= salines) and assign to new dataframe
c <- d %>%
  filter(Stimulant == "saline") %>%
  group_by(module_name) %>%
  summarise(mean_sal = mean(moduleTPM),
            sd_sal = sd(moduleTPM))

# calculate Z-scaled module scores
dat <- d %>%
  left_join(c) %>%
  mutate(module_Zscore = ((moduleTPM - mean_sal)/sd_sal))

### Correlation matrix (D7 samples)
# select D7 samples and reformat data
df <- dat %>%
  filter(Stimulant == "TST_D7") %>%
  filter(module_name %in% c("T_cell","myeloid","NK","CCND1","CD4_T","CD8_T")) %>%
  select(sample,module_name,module_Zscore) %>%
  pivot_wider(names_from = module_name, values_from = module_Zscore) %>%
  column_to_rownames("sample") %>%
  rename('All T cells' = T_cell,
         'CD4 T cells' = CD4_T,
         'CD8 T cells' = CD8_T,
         'NK cells' = NK,
         'Myeloid cells' = myeloid,
         'Proliferation module' = CCND1)

# generate correlation matrix
cor.mat <- as.matrix(round(cor(df, use = "pairwise.complete.obs", method = "spearman"),5))

# plot heatmap
col_fun = colorRamp2(c(0,0.5,1), c("#0270A0", "#E6EAC0", "#950018"))
cm = "average"

hm <- Heatmap(cor.mat,
              border_gp = gpar(col = "black", lty = 1),
              show_row_names = T,
              row_names_gp = gpar(fontsize = 10, fontface="bold"),
              show_column_names = F,
              col = col_fun,
              cluster_rows = T,
              cluster_columns = T,
              clustering_method_rows = cm,
              clustering_method_columns = cm,
              heatmap_legend_param = list(
                direction = "horizontal",
                at = c(0,0.5,1),
                labels = c(0,0.5,1),
                labels_gp = gpar(fontsize = 10, fontface="bold"),
                title = "Spearman correlation",
                legend_width = unit(2, "cm"),
                title_position = "lefttop",
                title_gp = gpar(fontsize = 10, fontface="bold"))
)
draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right")
