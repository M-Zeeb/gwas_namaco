############################################################
# Project: SHCGWAS
# Script:  heritability_mixed.R
# Author:  Marius Zeeb
# Date:    2025-11-10
# Purpose: Calculates heritability using linear mixed models
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

pheno <- args

library(tidyverse)
library(foreach)
library(doParallel)
library(ape)
library(phangorn)
library(lme4)

pheno_data <- read.csv("./phenotype_files/phenotype_file.csv")

heri_covariables <- readLines("./phenotype_files/phenotype_covariables.txt")

# ---------------------------
# 2. Functions
# ---------------------------

filter_tree <- function(lineage_tree, pheno, pheno_data) {
  ## Filter tree
  ## drop tips
  ## read pheno type data and select those present

  ## Take longest sequence (least gaps) per patient in case of multiple sequences
  ## (last resort, ideally done before)
  alignment_for_filter <- Biostrings::readAAStringSet(paste0("./alignment_nt/aligned.fa"))

  alignment_for_filter <- DECIPHER::RemoveGaps(alignment_for_filter)

  alignment_for_filter <- as.data.frame(cbind(
    alignment_for_filter@ranges@NAMES,
    as.numeric(alignment_for_filter@ranges@width)
  ))

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

  ## Drop sequences without pheno data
  lineage_tree <- drop.tip(
    lineage_tree,
    lineage_tree$tip.label[!(lineage_tree$tip.label %in% pheno_data$id |
      lineage_tree$tip.label == "Ref.D.TZ.01.A280.AY253311"
    )]
  )

  lineage_tree <- root(lineage_tree, "Ref.D.TZ.01.A280.AY253311")

  return(lineage_tree)
}

cluster_for_xtt <- function(subtres) {
  ## calculated maxdistance and number of tips for each subtree
  cluster_max_dist <- foreach(x = 1:length(subtres), .combine = "rbind") %do% {
    return(cbind(x, max(cophenetic.phylo(subtres[[x]])), subtres[[x]]$Ntip))
  }

  ## data formatting
  cluster_max_dist <- as.data.frame(cluster_max_dist)
  colnames(cluster_max_dist) <- c("subtree", "max_dist", "ntips")

  return(cluster_max_dist)
}

extract_dist_clusters <- function(cluster_max_dist, dist_thres, subtres) {
  ## filter to clusters with a within maximum distance of xx
  cluster_max_dist_for_filter <- cluster_max_dist[cluster_max_dist$max_dist < dist_thres, ]

  if (nrow(cluster_max_dist_for_filter) == 0) {
    return(NULL)
  }

  ## extract maximum cluster
  clusters <- NULL
  for (i in 1:nrow(cluster_max_dist_for_filter)) {
    # end loop when all clusters are checked
    if (nrow(cluster_max_dist_for_filter) == 0) {
      break
    }

    ## index of largest cluster (take first one if multiple with same size)
    largest_clus_index <- which(cluster_max_dist_for_filter$ntips == max(cluster_max_dist_for_filter$ntips))[1]

    ## clustername of largest cluster (index of subtree)
    cluster_name_maxtips_temp <- cluster_max_dist_for_filter$subtree[largest_clus_index]

    ## tip names in selected cluster
    tipsnames <- subtres[[cluster_name_maxtips_temp]]$tip.label

    ## remove cluster from (no longer needed)
    cluster_max_dist_for_filter <- cluster_max_dist_for_filter[-largest_clus_index, ]

    ## tips already accounted for?
    ## if yes, next start with next cluster
    ## possible can be more efficient with recursion to filter out all already accounted clusters directly
    ## but its fast enough I would say
    if (sum(tipsnames %in% clusters[, 2]) > 0) {
      next
    }

    ## bind cluster name and tipnames within cluster
    clusters <- rbind(
      clusters,
      cbind(cluster_name_maxtips_temp, tipsnames)
    )
  }

  clusters <- as.data.frame(clusters) %>%
    group_by(cluster_name_maxtips_temp) %>%
    mutate(n = n())

  return(clusters)
}

xtt_data_general <- function(pheno, pheno_data) {
  tree <- read.tree(paste0("./phylos/aligned.fa.treefile"))

  lineage_tree <- filter_tree(tree, pheno, pheno_data)

  lineage_tree <- drop.tip(
    lineage_tree,
    lineage_tree$tip.label[lineage_tree$tip.label == "Ref.D.TZ.01.A280.AY253311"]
  )

  ## extract all subtrees
  subtres <- subtrees(lineage_tree)

  tree_clusters <- cluster_for_xtt(subtres)

  for (dist_thres in as.numeric(c(seq(0.003, 0.009, 0.001), seq(0.01, 0.09, 0.01), seq(0.1, 0.5, 0.1)))) {
    clusters_extracted <- extract_dist_clusters(tree_clusters, dist_thres, subtres)
    if (is.null(clusters_extracted)) {
      next
    }

    clusters_dists <- clusters_extracted %>% select(id = tipsnames, cluster = cluster_name_maxtips_temp)

    clusters_dists <- clusters_dists %>%
      group_by(cluster) %>%
      mutate(n = n()) %>%
      filter(n > 1)

    clusters_dists <- clusters_dists %>%
      select(id, cluster)

    write.csv(
      clusters_dists,
      paste0(
        "./phylo_clusters/",
        pheno, "_", dist_thres, ".csv"
      )
    )
  }
}

heritability_lmer <- function(lmer_fit, alpha = 0.05) {
  # Extract variance components
  vc <- as.data.frame(VarCorr(lmer_fit))
  psi2 <- vc$vcov[vc$grp == "cluster"] # random-effect variance
  sigma2 <- attr(VarCorr(lmer_fit), "sc")^2 # residual variance

  h2 <- psi2 / (psi2 + sigma2)

  set.seed(123)

  # Bootstrap heritability estimate
  boot_fun <- function(fit) {
    vc <- as.data.frame(VarCorr(fit))
    psi2 <- vc$vcov[vc$grp == "cluster"]
    sigma2 <- attr(VarCorr(fit), "sc")^2
    psi2 / (psi2 + sigma2)
  }

  boot_res <- bootMer(lmer_fit, boot_fun, nsim = 1000)
  h_ci <- quantile(boot_res$t, probs = c(alpha / 2, 1 - alpha / 2))

  output <- array(c(h2, h_ci),
    dim = c(1, 3),
    dimnames = list(
      c("heritability"),
      c("Estimate", paste(formatC(c(alpha / 2, 1 - alpha / 2) * 100, digits = 1, format = "f"),
        "%",
        sep = ""
      ))
    )
  )
  return(output)
}


# ---------------------------
# 3. Calculate heritability
# ---------------------------

# Calculate clusters
xtt_data_general(pheno = pheno, pheno_data = pheno_data)

# Calculate heritabilities
registerDoParallel(cores = 6L)
heritabilities <- foreach(i = c(0.1, 0.2, 0.09, seq(0.03, 0.06, 0.01)), .errorhandling = "remove", .combine = "rbind") %dopar% {
  clusters_mm <- read.csv(paste0(
    "./phylo_clusters/",
    pheno, "_", i, ".csv"
  ))

  clusters_mm <- clusters_mm %>%
    mutate(cluster_new = cur_group_id(), .by = "cluster") %>%
    select(id, cluster = cluster_new)

  pheno_data_tmp <- merge(pheno_data, clusters_mm, by = "id", all = FALSE)

  pheno_data_tmp$cluster <- paste0("cl_", pheno_data_tmp$cluster)
  pheno_data_tmp$cluster <- factor(pheno_data_tmp$cluster, levels = sort(unique(pheno_data_tmp$cluster)))

  formula <- formula(paste0(pheno, " ~ ", paste(c("1", "(1|cluster)"), collapse = " + ")))
  fm1 <- lme4::lmer(formula, data = pheno_data_tmp)
  res1 <- heritability_lmer(fm1)
  h1 <- res1[1]
  h1_low <- res1[2]
  h1_up <- res1[3]
  p1 <- lmerTest::ranova(fm1)$`Pr(>Chisq)`[2]

  formula <- formula(paste0(pheno, " ~ ", paste(c(heri_covariables, "(1|cluster)"), collapse = " + ")))
  fm2 <- lme4::lmer(formula, data = pheno_data_tmp)
  res2 <- heritability_lmer(fm2)
  h2 <- res2[1]
  h2_low <- res2[2]
  h2_up <- res2[3]
  p2 <- lmerTest::ranova(fm2)$`Pr(>Chisq)`[2]

  heris_1 <- c(
    as.character(pheno), i,
    length(unique(pheno_data_tmp[, "cluster"])), "unadjusted", h1, h1_low, h1_up, p1, nrow(pheno_data_tmp)
  )

  heris_2 <- c(
    as.character(pheno), i,
    length(unique(pheno_data_tmp[, "cluster"])), "adjusted", h2, h2_low, h2_up, p2, nrow(pheno_data_tmp)
  )

  return(rbind(heris_1, heris_2))
}

print(heritabilities)
if (!is.null(heritabilities)) {
  colnames(heritabilities) <- c(
    "phenotype", "cluster_threshold", "n_clusters", "adjusted",
    "heritability", "h_low", "h_up", "pval", "n"
  )

  write.csv(heritabilities, paste0("./results/", pheno, "_mixed_heri.csv"), row.names = FALSE)
}
