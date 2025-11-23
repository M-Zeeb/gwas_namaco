############################################################
# Project: SHCGWAS
# Script:  prepare_poumm.R
# Author:  Marius Zeeb
# Date:    2025-11-05
# Purpose: Performs data preparation four poumm heritability analysis
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

print(args)

pheno <- args

library(ape)
library(tidyverse)

lineage_tree <- read.tree("./phylos/aligned.fa.treefile")

alignment_for_filter <- Biostrings::readAAStringSet(paste0("./alignment_nt/aligned.fa"))

pheno_data <- read.csv("./phenotype_files/phenotype_file.csv")

heri_covariables <- readLines("./phenotype_files/phenotype_covariables.txt")

# ---------------------------
# 2. Prepare data
# ---------------------------

##Filter tree
##drop tips
##read pheno type data and select those present

##Take longest sequence (least gaps) per patient in case of multiple sequences
##(last resort, ideally done before)

alignment_for_filter <- DECIPHER::RemoveGaps(alignment_for_filter)

alignment_for_filter <- as.data.frame(
    cbind(
        alignment_for_filter@ranges@NAMES,
        as.numeric(alignment_for_filter@ranges@width)
        )
    ) %>%
    mutate(V1 = gsub("/.+", "", V1))

pheno_data <- pheno_data %>% 
    filter(!is.na(get(pheno)))

pheno_data <- merge(pheno_data,
                    alignment_for_filter,
                    by.x = "id",
                    by.y = "V1",
                    all = FALSE
                    )

pheno_data <- pheno_data %>%
    mutate(V2 = as.numeric(V2)) %>%
    group_by(id) %>%
    arrange(desc(V2), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

##Drop sequences without pheno data
lineage_tree$tip.label = gsub("/.+", "", lineage_tree$tip.label)
lineage_tree <- drop.tip(lineage_tree,
                        lineage_tree$tip.label[!(lineage_tree$tip.label %in% pheno_data$id |
                                                lineage_tree$tip.label == "Ref.D.TZ.01.A280.AY253311"
                                                )
                                            ]
                        )

lineage_tree <- root(lineage_tree, "Ref.D.TZ.01.A280.AY253311")

#save phenotype data
pheno_data <- pheno_data %>%
    filter(id %in% lineage_tree$tip.label)


formula <- formula(paste0(pheno, " ~ ", paste(heri_covariables, collapse = " + ")))
lm_res <- lm(formula, data = pheno_data)
pheno_data[, paste0(pheno, "_res")] <- lm_res$residuals


pheno_data <- pheno_data %>%
    mutate(fid = 0) %>%
    select(fid, id, pheno, paste0(pheno, "_res"))

#save filtered tree
write.tree(lineage_tree, paste0("./phylos/tree_poumm_", pheno, ".treefile"))

write.csv(pheno_data,
          paste0("./phenotype_files/poumm_", pheno, ".phen"),
          row.names = FALSE
        )