############################################################
# Project: SHCGWAS
# Script:  gwas_visualization.R
# Author:  Marius Zeeb
# Date:    2025-11-13
# Purpose: Creates GWAS figure
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

library(tidytable)
library(stringr)
library(stringi)
library(forcats)
library(lubridate)
library(svMisc)
library(ggrepel)
library(qdapTools)
library(foreach)
library(doParallel)

registerDoParallel(cores = 8L)

# ---------------------------
# 2. Functions
# ---------------------------

## manhattan plot
manhattan_gg <- function(data, multiple_testing_threshold) {

  data <- data %>%
    mutate(
        estimate = ifelse(estimate >= 0, "increase", "decrease"),
        region = stri_trans_totitle(region),
        dna = as.numeric(dna)
    )

  data %>%
    mutate(p.value = -log10(p.value)) %>%
    mutate(size_p = (exp(p.value) / exp(max(p.value))) * 3) %>%
    mutate(alpha_p = (p.value / max(p.value)) * 1) %>%
    ggplot(.) +
        geom_point(aes(
        x = dna,
        y = p.value,
        size = size_p,
        colour = pheno,
        shape = region
        )) +
        xlab(paste0("HIV genome [Base position]")) +
        ylab("Significance [-log10(pvalue)]") +
        geom_hline(yintercept = multiple_testing_threshold) +
        scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 9, 10, 11)) +
        theme_minimal() +
        labs(colour = "", shape = "") +
        ylim(0, 6.1) +
        guides(
        size = FALSE,
        colour = guide_legend(override.aes = list(size = 5)),
        shape = guide_legend(override.aes = list(size = 5))
        ) +
        theme(text = element_text(size = 40)) +
        geom_text_repel(aes(
        x = dna,
        y = p.value,
        label = ifelse(p.value > 3.5, paste0(ref, hxb2_index, term), "")
        ), size = 6)
}

## manhattan plot visual abstract
manhattan_gg_visual <- function(data, multiple_testing_threshold) {
    data <- data %>%
    mutate(
        estimate = ifelse(estimate >= 0, "increase", "decrease"),
        region = stri_trans_totitle(region),
        dna = as.numeric(dna)
    )

  # data = merge(data_max,data_auc, by = c("p.value","type","hxb2_index_env","estimate"), all = TRUE)
  data %>%
    mutate(p.value = -log10(p.value)) %>%
    mutate(size_p = (exp(p.value) / exp(max(p.value))) * 3) %>%
    mutate(alpha_p = (p.value / max(p.value)) * 1) %>%
    ggplot(.) +
    geom_point(aes(
      x = dna,
      y = p.value,
      size = size_p,
      colour = pheno
    )) +
    xlab(paste0("HIV genome")) +
    ylab("Significance") +
    geom_hline(yintercept = multiple_testing_threshold) +
    theme_minimal() +
    labs(colour = "", shape = "") +
    ylim(0, 6.1) +
    guides(size = FALSE, colour = guide_legend(override.aes = list(
      size =
        5
    ))) +
    theme(
      axis.title.x = element_text(size = 35),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20),
      text = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.text = element_text(size = 20),
      panel.spacing = unit(2, "lines"),
      strip.text.x = element_text(size = 35)
    ) #+
  # geom_text_repel(aes(x = dna, y = p.value,label = ifelse(p.value>3.5,paste0(ref,hxb2_index,term),'')),size = 6)
}

## effective test size for helper function for "outer"-function
get.V2 <- function(x, y) {
  rcompanion::cramerV(x, y, bias.correct = TRUE)
}

## effective test size
effective_test_size <- function(subtype, data, file_annote) {
  # subtype = "B"
  #  data = gwas_res_B_clean
  #  region_indx = "gag"
  eff_tests <- foreach(
    region_indx = c("env", "tat", "gag", "pol", "rev", "vif", "vpr", "vpu", "nef"),
    .combine = "rbind"
  ) %dopar% {
    eff_tes_siz <- Biostrings::readAAStringSet(
      paste0(
        "/Users/mariuszeeb/Documents/neuro/alignment/",
        region_indx,
        "/nci_AA_",
        subtype,
        "_",
        region_indx,
        "_aligned_",
        file_annote,
        ".fa"
      )
    )

    for_cramer <- as.matrix(eff_tes_siz[-1, ])

    colnames(for_cramer) <- paste0("pos", seq(1:ncol(for_cramer)))

    for_cramer <- for_cramer[, colnames(for_cramer) %in% unique(paste0("pos", data$position[data$region == region_indx]))]

    cramers_v <- v_outer(for_cramer, get.V2)

    eigen <- as.data.frame(eigen(cramers_v)$values)
    colnames(eigen) <- "eigen"

    effec_test <- eigen %>%
      mutate(eigen = abs(eigen)) %>%
      dplyr::arrange(desc(eigen)) %>%
      mutate(n = dplyr::row_number()) %>%
      mutate(totaleigen = sum(eigen)) %>%
      mutate(eigenthreshold = totaleigen * 0.995) %>%
      mutate(eigensum = cumsum(eigen)) %>%
      mutate(eigensum_perc = eigensum / totaleigen) %>%
      mutate(effect_tests = eigensum >= eigenthreshold) %>%
      filter(effect_tests == TRUE) %>%
      select(n) %>%
      min()

    eff_test <- length(data$position[data$region == region_indx]) * (as.numeric(effec_test) /
      as.numeric(ncol(for_cramer)))

    return(eff_test)
  }

  return(eff_tests)
}

# ---------------------------
# 3. Visualize
# ---------------------------

which_gwas <- grep("gwas_clean", list.files("./results"), value = TRUE)
gwases <- lapply(which_gwas, function(g) read.csv(paste0("./results/", g)))
gwases <- bind_rows(gwases)

welp <- manhattan_gg(gwases, -log10(0.05))

ggsave(welp, filename = paste0("results/gwas_manhattan.png"), width = 30, height = 15)
