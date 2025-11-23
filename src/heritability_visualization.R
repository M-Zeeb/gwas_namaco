############################################################
# Project: SHCGWAS
# Script:  heritability_visualization.R
# Author:  Marius Zeeb
# Date:    2025-11-13
# Purpose: Creates heritability figure
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

library(tidyverse)
library(foreach)
library(doParallel)

phenos <- readLines("phenotype_files/phenotype_outcome.txt")

# ---------------------------
# 2. Mixed Heritability
# ---------------------------

heritabilities <- lapply(phenos, function(pheno) read.csv(paste0("./results/", pheno, "_mixed_heri.csv")))
heritabilities <- bind_rows(heritabilities)

heritabilities_forplot <- heritabilities %>%
    as.data.frame() %>%
    mutate(
        cluster_threshold = factor(cluster_threshold,
            levels = c("0.1", "0.09", "0.08", "0.07", "0.06", "0.05", "0.04", "0.03", "0.02", "0.01")
        ),
        pos = ifelse(adjusted == "adjusted", -0.01, 0.01),
        Adjusted = ifelse(adjusted == "adjusted", "Multivariable", "Univariable"),
        Adjusted = factor(Adjusted, levels = c("Univariable", "Multivariable")),
        h_low = as.numeric(h_low),
        h_up = as.numeric(h_up),
        heritability = as.numeric(heritability),
        h_low = h_low * 100,
        h_up = h_up * 100,
        heritability = heritability * 100
    )

heri_plot <- ggplot(heritabilities_forplot, aes(color = Adjusted)) +
    geom_errorbar(aes(xmin = h_low, xmax = h_up, y = pos), linewidth = 1.1, width = 0.0175) +
    geom_point(aes(x = heritability, y = pos), size = 4) +
    labs(
        color = "",
        y = "Phylogenetic cluster\ndistance threshold",
        x = "Broad Sense Heritability [%]"
    ) +
    theme_light() +
    theme(
        legend.position = "right", text = element_text(size = 20),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(size = 15)
    ) +
    scale_color_manual(values = c("#833ab4", "#fd1d1d")) +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "1"), limits = c(-20, 100)) +
    geom_vline(xintercept = 0) +
    facet_grid(cluster_threshold ~ phenotype)

ggsave(
    heri_plot,
    filename = "results/heritability_mixed_fig.png",
    width = 15,
    height = 4
)

# ---------------------------
# 3. POUMM Heritability
# ---------------------------

registerDoParallel(6L)

for (pheno in phenos) {
    print(pheno)

    which_results <- list.files(paste0("poumm/", pheno, "_bs"))

    foreach(i = which_results) %dopar% {
        results_raw <- readRDS(paste0("poumm/", pheno, "_bs/", i))

        results <- summary(results_raw)

        results$pheno <- pheno
        results$res <- ifelse(grepl("res", i), "res", NA)
        results$bs_index <- sub(".*_(\\d+)\\.RData$", "\\1", i)

        ci_heri <- unlist(results$HPD)
        dim(ci_heri) <- c(2, nrow(results))
        results$lower <- t(ci_heri)[, 1]
        results$upper <- t(ci_heri)[, 2]

        results <- results %>%
            as.data.frame() %>%
            select(-HPD)

        write.table(
            results,
            paste0(
                "poumm/bs_summaries_",
                ifelse(grepl("pmm", i), "pmm", "poumm"),
                "/",
                gsub("RData", "csv", i)
            ),
            row.names = FALSE,
            col.names = FALSE
        )
    }
}

system("cat poumm/bs_summaries_poumm/* > poumm/poumm_results.csv")
system("cat poumm/bs_summaries_pmm/* > poumm/pmm_results.csv")


bs_results_pmm <- read.table("poumm/pmm_results.csv")
bs_results_poumm <- read.table("poumm/poumm_results.csv")

colnames(bs_results_poumm) <- colnames(bs_results_pmm) <- c("stat", "N", "MLE", "PostMean", "ESS", "G.R.", "pheno", "res", "bs_index", "lower", "upper")

bs_results_pmm <- bs_results_pmm %>%
    filter(stat == "H2e") %>%
    mutate(res = ifelse(is.na(res), "univar", res)) %>%
    group_by(pheno, res) %>%
    mutate(
        lower_ci = mean(lower),
        upper_ci = mean(upper),
        mean = mean(PostMean),
        model = "pmm"
    ) %>%
    slice(1)

bs_results_poumm <- bs_results_poumm %>%
    filter(stat == "H2e") %>%
    mutate(res = ifelse(is.na(res), "univar", res)) %>%
    group_by(pheno, res) %>%
    mutate(
        lower_ci = mean(lower),
        upper_ci = mean(upper),
        mean = mean(PostMean),
        model = "poumm"
    ) %>%
    slice(1)

bs_results <- rbind(bs_results_pmm, bs_results_poumm) %>%
    mutate(
        pos_1 = ifelse(model == "poumm", 14, 4),
        pos_2 = ifelse(res == "res", -0.75, 0.75)
    ) %>%
    ungroup() %>%
    mutate(
        adjusted = ifelse(res == "res", "Multivariable", "Univariable"),
        adjusted = factor(adjusted, levels = c("Univariable", "Multivariable")),
        pos = pos_1 + pos_2
    )

heri_plot <- bs_results %>%
    ggplot(., aes(color = adjusted)) +
    geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci, y = pos_1 + pos_2), linewidth = 1.1, width = 0.5) +
    geom_point(aes(x = mean, y = pos_1 + pos_2), size = 4) +
    scale_y_continuous(
        breaks = c(
            4,
            14
        ),
        labels = c(
            "BM",
            "OU"
        )
    ) +
    labs(color = "", y = "", x = "Broad Sense Heritability") +
    theme_light() +
    theme(
        legend.position = "bottom", text = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12)
    ) +
    scale_color_manual(values = c("#833ab4", "#fd1d1d")) +
    scale_alpha_manual(values = c(0.5, 1)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"), limits = c(-1, 1)) +
    geom_vline(xintercept = 0) +
    facet_grid(~pheno)

ggsave(heri_plot, filename = paste0("results/heritability_poumm_fig.png"), width = 15, height = 5)
