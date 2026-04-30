############################################################
# Project: SHCGWAS
# Script:  gwas.R
# Author:  Marius Zeeb
# Date:    2025-09-26
# Purpose: Performs GWAS
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
library(brglm2)
registerDoParallel(cores = 8L)

# ---------------------------
# 2. Functions
# ---------------------------

lm_custom <- function(outcome, covariables, data_pos) {
  get_vars_to_omit <- \(d){
    names(d) |>
      sapply(FUN = \(name){
        xs <- d[[name]]
        is.numeric(xs) | (length(levels(as.factor(as.character(xs)))) > 1)
      }) |>
      Filter(f = \(xs) !xs) |>
      names()
  }

  covariables_take <- covariables[!covariables %in% get_vars_to_omit(data_pos)]
  formula_covars_update <- formula(paste0(
    outcome, " ~ fct_infreq(aa) + ",
    paste(covariables_take, collapse = "+")
  ))

  return(lm(formula_covars_update, data = data_pos))
}

glm_custom <- function(outcome, covariables, data_pos) {
  get_vars_to_omit <- \(d){
    names(d) |>
      sapply(FUN = \(name){
        xs <- d[[name]]
        is.numeric(xs) | (length(levels(as.factor(as.character(xs)))) > 1)
      }) |>
      Filter(f = \(xs) !xs) |>
      names()
  }

  covariables_take <- covariables[!covariables %in% get_vars_to_omit(data_pos)]
  formula_covars_update <- formula(paste0(
    outcome, " ~ fct_infreq(aa) + ",
    paste(covariables_take, collapse = "+")
  ))

  return(glm(formula_covars_update, data = data_pos, family = "binomial", method = brglmFit))
}

gwas_lm_multi_var <- function(gwas_single_test_set, outcomes, covariables) {
  result_all_outcomes <- foreach(outcome = outcomes, .combine = "rbind") %do% {
    result_outcome <- gwas_single_test_set %>%
      dplyr::group_by(position) %>%
      dplyr::do(lm = broom::tidy(
        lm_custom(
          outcome,
          covariables,
          data_pos = .
        )
      )) %>%
      unnest(lm) %>%
      mutate(pheno = outcome)

    return(result_outcome)
  }

  return(result_all_outcomes)
}

gwas_glm_multi_var <- function(gwas_single_test_set, outcomes, covariables) {
  result_all_outcomes <- foreach(outcome = outcomes, .combine = "rbind") %do% {
    result_outcome <- gwas_single_test_set %>%
      filter(sum(get(outcome)) > 1, .by = c(position, aa)) %>%
      filter(length(unique(aa)) > 1, .by = position) %>%
      dplyr::group_by(position) %>%
      dplyr::do(glm = broom::tidy(
        glm_custom(
          outcome,
          covariables,
          data_pos = .
        )
      )) %>%
      unnest(glm) %>%
      mutate(pheno = outcome)

    return(result_outcome)
  }

  return(result_all_outcomes)
}

filter_aa_results <- function(results, raw_data, seqs_to_keep) {
  results <- results %>%
    filter(grepl("aa", term)) %>%
    mutate(position = gsub("pos", "", position)) %>%
    mutate(term = gsub("fct\\_infreq\\(aa\\)", "", term)) %>%
    mutate(interaction = grepl(":gag", term)) %>%
    mutate(term = gsub(":gag", "", term))

  results <- results %>%
    mutate(position = as.numeric(position)) %>%
    full_join(.,
      seqs_to_keep[, c("msa_pos", "hxb2_index")],
      by = c("position" = "msa_pos")
    ) %>%
    arrange(pheno, position) %>%
    fill(hxb2_index, .direction = "down", .by = pheno)

  results <- results[complete.cases(results), ]

  results <- raw_data %>%
    mutate(position = gsub("pos", "", position)) %>%
    mutate(position = as.numeric(position)) %>%
    slice(1, .by = c("position", "aa")) %>%
    select(position, aa, n, ref, n_ref) %>%
    left_join(
      results,
      .,
      by = c("position" = "position", "term" = "aa")
    )

  return(results)
}

# ---------------------------
# 3. RUN GWAS
# ---------------------------

gwas_AA_set <- read.csv(paste0("gwas_data/", region, "_gwas_data.csv"))
gwas_AA_idx <- read.csv(paste0("gwas_data/", region, "_seq_idx.csv"))

gwas_outcomes <- readLines("phenotype_files/phenotype_outcome.txt")
gwas_covariables <- readLines("phenotype_files/phenotype_covariables.txt")#[c(1,2,15,16,17)] # for reduced covariables

gwas_res <- foreach(outcome = gwas_outcomes, .combine = "rbind", .errorhandling = "remove") %dopar% {

  if (length(unique(gwas_AA_set[, outcome])) == 2) {
    gwas_res_tmp <- gwas_glm_multi_var(gwas_AA_set, outcome, gwas_covariables)
  } else {
    gwas_res_tmp <- gwas_lm_multi_var(gwas_AA_set, outcome, gwas_covariables)
  }

  if (nrow(gwas_res_tmp) == 0) {
    warning(paste0("No results for outcome: ", outcome))
    return(NULL)
  }
  return(gwas_res_tmp)
}


# ---------------------------
# 4. Clean GWAS
# ---------------------------

## clean gwas results
position_hxb2_ref <- list(
  # env
  c(seq(6046, 8794, 3)),
  # gag
  c(seq(791, 2291, 3)),
  # pol
  c(seq(2066, 5095, 3)),
  # tat
  c(c(seq(5831, 6045, 3), seq(8379, 8469, 3))),
  # nef
  c(seq(8797, 9417, 3)),
  # rev
  c(c(seq(5970, 6045, 3), seq(8381, 8653, 3))),
  # vpu
  c(seq(6062, 6310, 3)),
  # vif
  c(seq(5041, 5619, 3)),
  # vpr
  c(c(seq(5559, 5771, 3), seq(5773, 5850, 3))),
  # whole
  c(seq(1, 9719))
)
names(position_hxb2_ref) <- c("env", "gag", "pol", "tat", "nef", "rev", "vpu", "vif", "vpr", "whole")

posis <- as.data.frame(cbind(position_hxb2_ref[[region]], 1:length(position_hxb2_ref[[region]])))
colnames(posis) <- c("dna", "aa")

gwas_res_clean <- filter_aa_results(gwas_res, gwas_AA_set, gwas_AA_idx) %>%
  inner_join(
    .,
    posis,
    by = c("hxb2_index" = "aa")
  ) %>%
  mutate(region = region)

write.csv(gwas_res_clean, paste0("results/", region, "_gwas_clean.csv"), row.names = FALSE)

print(sessionInfo())