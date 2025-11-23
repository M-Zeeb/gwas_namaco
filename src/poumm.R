############################################################
# Project: SHCGWAS
# Script:  poumm.R
# Author:  Marius Zeeb
# Date:    2025-11-05
# Purpose: Performs Heritability analysis with poumm with bootstrap
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

print(args)

pheno <- args

library(POUMM)
library(ape)
library(phylolm)
library(phangorn)
library(tidyverse)
library(doParallel)

registerDoParallel(cores = 6)

# create folder
ifelse(!dir.exists(file.path(paste0("./poumm/",pheno,"_bs"))), dir.create(file.path(paste0("./poumm/",pheno,"_bs"))), FALSE)

# read ultrametric tree
lineage_tree_bs <- read.tree(paste0("./phylos/aligned.fa.ufboot"))


# ---------------------------
# 2. RUN POUMM
# ---------------------------

for(i in 1:3){

    lineage_tree <- lineage_tree_bs[[i]]

    ##read phenotype
    pheno_data <- as.data.frame(read.csv(paste0("./phenotype_files/poumm_", pheno, ".phen"), header = TRUE))

    colnames(pheno_data) <- c("fid", "id", "pheno", "pheno_res")

    ##drop tips again to be safe
    lineage_tree <- drop.tip(
        lineage_tree,
        lineage_tree$tip.label[!lineage_tree$tip.label %in% pheno_data$id]
    )

    ##filter pheno to make sure to include only the ones in tree
    pheno_data <- pheno_data %>%
        filter(id %in% lineage_tree$tip.label) %>%
        arrange(factor(id, levels = lineage_tree$tip.label))

    lineage_tree <- reorder.phylo(lineage_tree, "pruningwise")
    lineage_tree_bm <- reorder.phylo(lineage_tree, "pruningwise")

    n <- length(lineage_tree$tip.label)
    D <- numeric(n)
    des <- lineage_tree$edge[, 2]
    externalEdge <- (des <= n)
    root.to.tip <- phylolm::pruningwise.distFromRoot(lineage_tree)[1:n]

    if (!is.ultrametric(lineage_tree)) {
        flag <- 1
        dis <- phylolm::pruningwise.distFromRoot(lineage_tree) # all root to tip distances, and root to node distances
        names(dis) <- lineage_tree$tip.label#??number of tip labels less than number of dis
        D <- max(dis[1:n]) - dis[1:n]
        D <- D - mean(D) #0-mean centered
        lineage_tree$edge.length[externalEdge] <- lineage_tree$edge.length[externalEdge] + D[des[externalEdge]]
    }

    ##Fix negative edge lengths
    lineage_tree <- nnls.tree(cophenetic(lineage_tree), lineage_tree, rooted=TRUE)

    ##Fix zero edge lengths
    lineage_tree$edge.length <- lineage_tree$edge.length + 0.000001

    #pheno$pheno <- (pheno$pheno - mean(pheno$pheno))/sd(pheno$pheno)

    ##phylogenetic mixed model comparision (is poumm better)? (take poumm_res from main poumm model below)
    #univariable
    specPMM <- specifyPMM(
        pheno_data$pheno,
        lineage_tree_bm,
        nSamplesMCMC = 9e5,
        parallelMCMC = TRUE
    )
    fitPMM <- POUMM(
        pheno_data$pheno,
        lineage_tree_bm,
        spec = specPMM,
        doMCMC = TRUE
    )

    specH2tMean <- specifyPOUMM_ATH2tMeanSeG0(
        pheno_data$pheno,
        lineage_tree,
        nSamplesMCMC = 9e5,
        parallelMCMC = TRUE
    )
    poumm_res <- POUMM::POUMM(
        pheno_data$pheno,
        lineage_tree, spec = list(nSamplesMCMC = 9e5, parallelMCMC = TRUE)
    )

    #Multivariable
    specPMM_multi <- specifyPMM(
        pheno_data$pheno_res,
        lineage_tree_bm,
        nSamplesMCMC = 9e5,
        parallelMCMC = TRUE
    )
    fitPMM_multi <- POUMM(
        pheno_data$pheno_res,
        lineage_tree_bm,
        spec = specPMM_multi,
        doMCMC = TRUE
    )

    specH2tMean_multi <- specifyPOUMM_ATH2tMeanSeG0(
        pheno_data$pheno_res,
        lineage_tree,
        nSamplesMCMC = 9e5,
        parallelMCMC = TRUE
    )
    poumm_res_multi <- POUMM::POUMM(
        pheno_data$pheno_res,
        lineage_tree,
        spec = list(nSamplesMCMC = 9e5, parallelMCMC = TRUE)
    )


    pheno_data$id <- as.character(pheno_data$id)

    saveRDS(poumm_res, file = paste0("./poumm/", pheno, "_bs/poumm_", pheno, "_", i, ".RData"))
    saveRDS(fitPMM, file = paste0("./poumm/", pheno, "_bs/pmm_", pheno, "_", i, ".RData"))
    saveRDS(poumm_res_multi, file = paste0("./poumm/", pheno, "_bs/poumm_", pheno, "_res_", i, ".RData"))
    saveRDS(fitPMM_multi, file = paste0("./poumm/", pheno, "_bs/pmm_", pheno, "_res_", i, ".RData"))

    print(paste0(i, " done"))
}