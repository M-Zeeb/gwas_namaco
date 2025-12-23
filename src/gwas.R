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
registerDoParallel(cores = 8L)

# ---------------------------
# 2. Functions
# ---------------------------

gwas_logistic_multi_var <- function(gwas_single_test_set) {
  # gwas_single_test_set = all_gwas[["env"]][[1]]

  ## filter out "dis" with single factor level
  gwas_single_test_set <- gwas_single_test_set %>%
    group_by(position) %>%
    na.omit() %>%
    mutate(n_dis = length(unique(dis))) %>%
    filter(n_dis > 1)
  # #
  # for(i in unique(gwas_single_test_set$position)){
  #   print(i)
  #
  #   hmpf = gwas_single_test_set %>% filter(position == i) %>% na.omit() %>%
  #     AER::tobit(auc_all_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .) %>% summary()
  #
  # }

  all_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(
      glm(
        auc_all_avg > 0 ~ fct_infreq(aa) + sex + age + education + risk + auc_rna +
          auc_cd4 + ethnicity + AUC_efv + antidepress + depri + drugs + dis + hep_c +
          hep_b + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
        data = .
      )
    )) %>%
    unnest(lm) %>%
    mutate(pheno = "all")

  slow_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(
      glm(
        auc_slow_avg > 0 ~ fct_infreq(aa) + sex + age + education + risk + auc_rna +
          auc_cd4 + ethnicity + AUC_efv + antidepress + depri + drugs + dis + hep_c +
          hep_b + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
        data = .
      )
    )) %>%
    unnest(lm) %>%
    mutate(pheno = "slow")

  conc_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(
      glm(
        auc_conc_avg > 0 ~ fct_infreq(aa) + sex + age + education + risk + auc_rna +
          auc_cd4 + ethnicity + AUC_efv + antidepress + depri + drugs + dis + hep_c +
          hep_b + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
        data = .
      )
    )) %>%
    unnest(lm) %>%
    mutate(pheno = "conc")

  freq_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(
      glm(
        auc_freq_avg > 0 ~ fct_infreq(aa) + sex + age + education + risk + auc_rna +
          auc_cd4 + ethnicity + AUC_efv + antidepress + depri + drugs + dis + hep_c +
          hep_b + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
        data = .
      )
    )) %>%
    unnest(lm) %>%
    mutate(pheno = "freq")

  all_nci <- rbind(all_nci, slow_nci, conc_nci, freq_nci)

  return(all_nci)
}

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

  return(glm(formula_covars_update, data = data_pos, family = "binomial"))
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
      group_by(position, aa) %>%
      filter(sum(get(outcome)) > 1) %>%
      group_by(position) %>%
      filter(length(unique(aa)) > 1) %>%
      do(glm = broom::tidy(
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

sgwas_multi_var_unadjusted <- function(gwas_single_test_set) {
  # gwas_single_test_set = all_gwas[["env"]][[1]]

  ## filter out "dis" with single factor level
  gwas_single_test_set <- gwas_single_test_set %>%
    group_by(position) %>%
    na.omit() %>%
    mutate(n_dis = length(unique(dis))) %>%
    filter(n_dis > 1)
  # #
  # for(i in unique(gwas_single_test_set$position)){
  #   print(i)
  #
  #   hmpf = gwas_single_test_set %>% filter(position == i) %>% na.omit() %>%
  #     AER::tobit(auc_all_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .) %>% summary()
  #
  # }

  all_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(summary(
      AER::tobit(
        auc_all_avg ~ fct_infreq(aa) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
          PC8 + PC9 + PC10,
        data = .
      )
    )$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "all")

  slow_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(summary(
      AER::tobit(
        auc_slow_avg ~ fct_infreq(aa) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
          PC8 + PC9 + PC10,
        data = .
      )
    )$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "slow")

  conc_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(summary(
      AER::tobit(
        auc_conc_avg ~ fct_infreq(aa) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
          PC8 + PC9 + PC10,
        data = .
      )
    )$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "conc")

  freq_nci <- gwas_single_test_set %>%
    # mutate(auc_avg = (auc_avg + 0.4279322)*-1) %>%
    group_by(position) %>%
    do(lm = broom::tidy(summary(
      AER::tobit(
        auc_freq_avg ~ fct_infreq(aa) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
          PC8 + PC9 + PC10,
        data = .
      )
    )$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "freq")

  all_nci <- rbind(all_nci, slow_nci, conc_nci, freq_nci)

  return(all_nci)
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
gwas_covariables <- readLines("phenotype_files/phenotype_covariables.txt")#[c(1,2,15,16,17)]

gwas_res <- foreach(outcome = gwas_outcomes, .combine = "rbind") %dopar% {
  gwas_lm_multi_var(gwas_AA_set, outcome, gwas_covariables)
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
  c(c(seq(5559, 5771, 3), seq(5773, 5850, 3)))
)
names(position_hxb2_ref) <- c("env", "gag", "pol", "tat", "nef", "rev", "vpu", "vif", "vpr")

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
