############################################################
# Project: SHCGWAS
# Script:  create_example_data.R
# Author:  Marius Zeeb
# Date:    2025-12-22
# Purpose: Creates example data for the GWAS/Heritability pipelines
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

library(tidytable)

# ---------------------------
# 2. Create example frequency files (based on random HIV reference sequences)
# ---------------------------

seqs <- Biostrings::readDNAStringSet("utility/references_seqs/HIV1_REF_2023_whole_DNA.fasta")

for(i in sample(1:555, 50, replace = FALSE)){
  
  seq <- as.matrix(seqs)[i,]
  seq <- seq[seq != "-"]
  seq_df <- cbind(seq,1:length(seq))
  seq_df <- cbind(seq_df, 
                  sample(c(1000:2000), nrow(seq_df), replace = TRUE),
                  sample(seq(from = 0.01, to = 0.4, by = 0.001), nrow(seq_df), replace = TRUE),
                  sample(c("A", "G", "C", "T"), nrow(seq_df), replace = TRUE)
                  )
  seq_df <- seq_df %>% as.data.frame() %>%
    mutate( V4 = as.numeric(V4) * 100,
            V2 = as.numeric(V2),
            V3 = as.numeric(V3),
            V4 = replace(V4, V5 == seq, NA)) %>%
    select(POS = V2, REF = seq, ALT = V5, AF = V4, COV = V3) 
  
  
  write.csv(seq_df, paste0("input/raw_seqs/", ids::uuid(1)))
}

# ---------------------------
# 3. Create files containing name of outcome and covariables
# ---------------------------

writeLines(c("Z_dummy1","Z_dummy2"), 
           paste0("input/phenotype_files/phenotype_outcome.txt")
           )

writeLines(c("age","sex"), 
           paste0("input/phenotype_files/phenotype_covariables.txt")
)

# ---------------------------
# 4. Generate random values for outcome and covariables
# ---------------------------

uuids <- list.files(paste0("input/raw_seqs/"))

phenotype_file <- data.frame(
  "ID" = sample(10000:99999, length(uuids)), 
  "Z_dummy1" = rnorm(length(uuids)), 
  "Z_dummy2" = rnorm(length(uuids)), 
  "age" = rnorm(length(uuids), mean = 53, sd = 10),
  "sex" = ifelse(rbinom(50,1,0.2) == 1, "Female", "Male"),
  "base_uuid" = uuids
)

write.csv(phenotype_file, 
          paste0("input/phenotype_files/phenotype_file.csv"),
          row.names = FALSE)

write.table(phenotype_file %>% select(base_uuid, ID), 
          paste0("input/phenotype_files/pheno.names.txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE)
