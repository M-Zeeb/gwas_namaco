############################################################
# Project: SHCGWAS
# Script:  dummy_phenotypes.R
# Author:  Marius Zeeb
# Date:    2025-09-26
# Purpose: To briefly generate a outcome file with covariates for GWAS testing
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

library(tidytable)

# ---------------------------
# 2. RUN
# ---------------------------

seq_uuids <- list.files("/Users/mariuszeeb/Documents/projects/SHCGWAS/raw_seqs_namaco/") %>% gsub(".csv","",.)
namaco_meta <- read.csv("~/Documents/projects/neuro_namaco_gwas/sequences/NGS_meta_information/NGS_samples_list_all runs_NGS_R1R2_v458_071024.csv")
namaco_meta <- namaco_meta %>% 
  filter(base_uuid %in% seq_uuids) %>%
  select(base_uuid,SHCS.ID)
shcs_pat <- read.csv("~/Documents/SHCS/csv_2509/pat.csv")
shcs_pat <- shcs_pat %>% inner_join(., namaco_meta, by = c("ID" = "SHCS.ID")) %>%
  select(base_uuid,ID,BORN,ETHNICITY,RISK)
write.csv(shcs_pat,"/Users/mariuszeeb/Documents/projects/SHCGWAS/pheno_namaco/phenotype.csv", row.names = FALSE)

print(sessionInfo())