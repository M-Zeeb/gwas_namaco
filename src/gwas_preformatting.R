############################################################
# Project: SHCGWAS
# Script:  gwas_preformatting.R
# Author:  Marius Zeeb
# Date:    2025-09-26
# Purpose: Performs necessary preformatting steps for GWAS
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

region <- args[1]  

library(tidytable)
library(stringr)
library(forcats)
library(AER)
library(broom)
library(lubridate)
library(svMisc)
library(ape)
library(phangorn)
library(ggrepel)
library(phytools)
library(qdapTools)
library(foreach)
library(doParallel)
registerDoParallel(cores = 8L)

# ---------------------------
# 2. Functions
# ---------------------------

reformat_seq_to_csv <- function(msa){
  
  msa_df <- as.data.frame(msa)
  msa_df$V1 <- msa@ranges@NAMES
  ##Rename column names
  names(msa_df) <- c("seq","base_uuid")
  msa_df <- msa_df %>% 
    dplyr::select(base_uuid,seq)
  ##Dataset for sequencing data
  msa <- msa_df %>% 
    separate(seq, into = paste("V", 1:(max(nchar(msa_df$seq))), sep = ""), sep = "")
  
  return(msa)
}

ref_seq_pos <- function(msa,region){
  
  hxb2_region <- paste0("HXB2_",toupper(region))
  
  seqs_to_keep <- as.matrix(msa %>%
                              filter(base_uuid == hxb2_region))
  
  seqs_to_keep <- t(seqs_to_keep)[-1,,drop = FALSE] %>% as.data.frame()
  colnames(seqs_to_keep) <- "variant"
  seqs_to_keep$gwas_pos <- as.numeric(gsub("V", "", rownames(seqs_to_keep)))  
  
  seqs_to_keep  <- seqs_to_keep  %>%
    filter(variant != "-") %>%
    mutate(hxb2_index = 1:n()) %>%
    mutate(hxb2_index_overall = as.numeric(hxb2_index))
  
  seqs_to_keep$msa_pos <- seqs_to_keep$gwas_pos 
  
  return(seqs_to_keep)
}

remove_hxb2_seq <- function(msa, region) {
  
  hxb2_region <- paste0("HXB2_", toupper(region))
  
  #remove hxb2
  msa <- msa %>%
    filter(.[[1]] != hxb2_region)
  
  colnames(msa)[2:(ncol(msa))] = paste0("pos", seq(1, (ncol(msa) - 1)))
  
  return(msa)
}

##build combined covariable/pheno data frame
pheno_covar_data_formatting <- function(region,
                                       pheno_covar_data,
                                       msa,
                                       wga_pca = FALSE) {
  
  if(wga_pca) {
    region <- "whole"
  }
  
  pca <- read.csv(
    paste0("pca/pca_",region,"/",region,".evec"),
    sep = "",
    header = FALSE
  )[-1, ]
  
  colnames(pca) <- c(
    "base_uuid",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
    "famid"
  )
  
  data <- pheno_covar_data %>%
    inner_join(.,
          pca[, -12],
          by = c("base_uuid" = "base_uuid")) %>%
    filter(base_uuid %in% unlist(msa[,1]))

  return(data)
}

##take a single sequence for each patients (ideally done before, but as a last resort here)
##Sequence with least missing is prefered
one_pat_one_seq <- function(msa,covar_data){
  
  msa_filter <- inner_join(msa,covar_data[,c("ID","base_uuid")], by = c("base_uuid" = "base_uuid"))
  msa_filter$missings <- rowSums(msa_filter[,-1] == "Z")
  msa_filter <- msa_filter %>% 
    arrange(ID, missings) %>% 
    slice(1, .by = ID) %>%
    select(-ID,-missings)
  
  return(msa_filter)
  
}

##dataset region generation
region_set <- function(msa, covar_data, maf, minref) {
  
  gwas_single_test <- msa %>%
    pivot_longer(!base_uuid, names_to = "position", values_to = "aa")
  
  gwas_single_test <- inner_join(
    gwas_single_test,
    covar_data,
    by = c("base_uuid" = "base_uuid")
  )

  gwas_single_test <- gwas_single_test %>%
    filter(aa != "Z")
  
  mincount <- gwas_single_test %>% 
    dplyr::select(position, aa) %>% 
    tidyr::gather(position, aa) %>% 
    dplyr::count(position, aa, .drop = FALSE) %>% 
    filter(n >= maf)
  
  refcount <- mincount %>% 
    group_by(position) %>% 
    arrange(position, desc(n)) %>% 
    slice(1, .by = position) %>% 
    select(position, ref = aa, n_ref = n)
  
  pos_to_take <- gwas_single_test %>% 
    dplyr::select(position, aa) %>% 
    tidyr::gather(position, aa) %>% 
    dplyr::count(position, aa, .drop = FALSE) %>% 
    filter(n >= maf) %>% 
    mutate(n_aas = sum(n > 1), .by = position) %>% 
    filter(n_aas > 1)
  
  gwas_single_test <- gwas_single_test %>% 
    filter(paste0(position, aa) %in% paste0(pos_to_take$position, pos_to_take$aa))
  
  gwas_single_test <- inner_join(gwas_single_test,
                           mincount,
                           by = c("position" = "position", "aa" = "aa")) %>% 
    inner_join(., refcount, by = c("position" = "position")) %>% 
    filter(n_ref > minref)
  
  return(gwas_single_test)
}

##function to combined subfuctions
gwas_ready_formatting_final <- function(region,
                                       maf = 5,
                                       minref = 10,
                                       wga_pca = TRUE) {
  
  ##geno type formatting (genotype file and hxb2 position relation as two output tables)
  msa <- Biostrings::readAAStringSet(paste0("alignment_aa/alignment_",region,".fa"))

  msa <- reformat_seq_to_csv(msa)
  
  seq_to_keep <- ref_seq_pos(msa, region)
  
  msa <- remove_hxb2_seq(msa, region)
  
  #Replace gaps "-" with "Z"
  msa <- msa %>%
    mutate(across(2:ncol(msa), ~ data.table::fifelse(.x == "-", "Z", .x)))
  
  pheno_covar_data <- read.csv("phenotype_files/phenotype_file.csv")
  
  ##combined dataset pheno covar
  covar_data <- pheno_covar_data_formatting(region, pheno_covar_data, msa, wga_pca = TRUE)
  
  ##take a single sequence for each patients (ideally done before, but as a last resort here)
  ##Sequence with least missing is prefered
  msa <- one_pat_one_seq(msa,covar_data)
  covar_data <- covar_data %>% 
    filter(base_uuid %in% msa$base_uuid)
  ##combined dataset geno pheno covar
  gwas_final <- region_set(msa, covar_data, maf, minref)
  
  write.csv(gwas_final, paste0("gwas_data/",region,"_gwas_data.csv"), row.names = FALSE)
  write.csv(seq_to_keep, paste0("gwas_data/",region,"_seq_idx.csv"), row.names = FALSE)
  
}

# ---------------------------
# 3. RUN GWAS preformatting
# ---------------------------

gwas_ready_formatting_final(region,  maf = 5, minref = 10, wga_pca = TRUE)

print(sessionInfo())