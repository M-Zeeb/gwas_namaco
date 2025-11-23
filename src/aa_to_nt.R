############################################################
# Project: SHCGWAS
# Script:  aa_to_nt.R
# Author:  Marius Zeeb
# Date:    2025-09-22
# Purpose: Generates NT alignment from AA alignment
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

region <- args[1]      

library(tidytable)
library(svMisc)
library(foreach)
library(doParallel)
library(Rcpp)
library(qdapTools)
library(glmnet)
library(Biostrings, include.only = c("xscat", "writeXStringSet", "readDNAStringSet","readAAStringSet"))
library(fastDummies)

# ---------------------------
# 2. Define function
# ---------------------------

AAtoNT = function(region,
                  depth_threshold) {
  
  outputpath <- paste0("alignment_aa/")
  
  file_name <- paste0(paste(
    c("alignment", region),
    collapse = "_"
  ), ".fa")
  
  aa <- Biostrings::readAAStringSet(paste0(outputpath, "/", file_name))
  if (region == "vpr") {
    aa <- aa[!aa@ranges@NAMES == "HXB2_VPR"]
    aa@ranges@NAMES[aa@ranges@NAMES == "NONHXB2_VPR"] <- "HXB2_VPR"
  }
  
  if (region %in% c("rev", "tat")) {
    blastmeta_1 <- read.csv(paste0("blastmeta/", paste(
      c("blastmeta", region, "exon1"), collapse = "_"
    ), ".csv")) %>% na.omit()
    
    blastmeta_2 <- read.csv(paste0("blastmeta/", paste(
      c("blastmeta", region, "exon2"), collapse = "_"
    ), ".csv")) %>% na.omit()
  }
  
  if (!region %in% c("rev", "tat")) {
    blastmeta <- read.csv(paste0("blastmeta/", paste(
      c("blastmeta", region), collapse = "_"
    ), ".csv")) %>% na.omit()
  }
  
  nt_al_list <- list()
  seq_ids <- c()
  freq_aligned_all <- NULL
  
  for (uuid in 1:length(aa)) {

    uuid_name <- aa@ranges@NAMES[uuid]
    
    aa_df <- as.data.frame(aa[uuid])
    
    aa_tab <- aa_df %>%
      separate(x, into = paste("V", 1:(max(nchar(aa_df$x))), sep = ""), sep = "")
    
    aa_tab <- as.data.frame(t(aa_tab))
    colnames(aa_tab) <- uuid_name
    
    aa_tab$row_num_aa <- seq(1:nrow(aa_tab))
    
    aa_tab <- aa_tab %>% 
      dplyr::filter(get(uuid_name) != "-")
    
    
    ##NT seqs

    if (uuid_name != paste0("HXB2_", toupper(region))) {
      print(uuid_name)
      NT <- Biostrings::readBStringSet(paste0(
        "codon_align/",
        paste(c("codon_align", region,"NT/"), collapse = "_"),
        uuid_name,
        ".fa"
      ))
      NT <- NT[!grepl(region, NT@ranges@NAMES, ignore.case = TRUE)]
      NT_df <- as.data.frame(NT[1])
      
      #base frequency file
      NT_freq <- read.csv(paste0(
        "freq_files/",
        uuid_name,
        ".csv"
      ))
      NT_freq_df <- as.data.frame(NT_freq)

      NT_freq_df$uuid <- uuid_name
      
      NT_freq_cov <- NT_freq_df %>% 
        filter(depth > depth_threshold & (!is.na(depth))) %>%
        mutate(freq_indx = cur_group_id(), .by = position)
      
      ##consesus whole sequence
      fasta_raw <- as.matrix(readDNAStringSet(
        paste0(
          "fasta_files/",
          uuid_name,
          ".fa"
        )
      ))
      fasta_raw <- as.data.frame(t(fasta_raw))
      fasta_raw$fasta_indx <- 1:nrow(fasta_raw)
      fasta_raw_temp <- fasta_raw[!fasta_raw[, 1] == "N", ]
      fasta_raw_temp$freq_indx <- 1:nrow(fasta_raw_temp)
      fasta_raw <- merge(fasta_raw,
                        fasta_raw_temp[, c("fasta_indx", "freq_indx")],
                        by = "fasta_indx",
                        all = TRUE)
      
      ##exon managing for rev and tat
      if (region %in% c("rev", "tat")) {
        ##extract blasted positions
        ##exon1
        fasta_raw_1 = NULL
        if (uuid_name %in% blastmeta_1$QueryID) {
          start_blast_1 = blastmeta_1$Q.start[blastmeta_1$QueryID == uuid_name]
          stop_blast_1 = blastmeta_1$Q.end[blastmeta_1$QueryID == uuid_name]
          fasta_raw_1 = fasta_raw %>% filter(fasta_indx >= start_blast_1 &
                                               fasta_indx <= stop_blast_1)
          #fasta_raw_1 = as.matrix(cbind(c(start_blast_1:stop_blast_1),fasta_raw[c(start_blast_1:stop_blast_1)]))
        }
        ##exon2
        fasta_raw_2 = NULL
        if (uuid_name %in% blastmeta_2$QueryID) {
          start_blast_2 <- blastmeta_2$Q.start[blastmeta_2$QueryID == uuid_name]
          stop_blast_2 <- blastmeta_2$Q.end[blastmeta_2$QueryID == uuid_name]
          fasta_raw_2 <- fasta_raw %>% filter(fasta_indx >= start_blast_2 &
                                               fasta_indx <= stop_blast_2)
          #fasta_raw_2 = as.matrix(cbind(c(start_blast_2:stop_blast_2),fasta_raw[,c(start_blast_2:stop_blast_2)]))
        }
        
        fasta_raw <- rbind(fasta_raw_1, fasta_raw_2)
        
      }
      
      
      ##
      ##start/stop index of blasted sequence
      if (!region %in% c("rev", "tat")) {
        start_blast <- blastmeta$Q.start[blastmeta$QueryID == uuid_name]
        stop_blast <- blastmeta$Q.end[blastmeta$QueryID == uuid_name]
        
        ##extract sequence from original fasta
        fasta_raw <- fasta_raw %>% 
          filter(fasta_indx >= start_blast &
                 fasta_indx <= stop_blast)
        #fasta_raw = as.matrix(cbind(start_blast:stop_blast,fasta_raw[,c(start_blast:stop_blast)]))
      }
    }
    
    ##hxb2 reference
    if (uuid_name == paste0("HXB2_", toupper(region))) {
      #random sequence to get an hxb2 nucleotide codon aligned sequence, from the next one or otherwise the one before
      NT <- Biostrings::readAAStringSet(
        paste0(
          "references_seqs/hxb2_",
          region,
          ".fa"
        )
      )
      
      ##for vpr due to frame shift, maybe hxb2 reference is not necessary at all, but lets remove it here and take only the non-frameshift complete vpr sequence
      if (region == "vpr") {
        NT <- NT[!NT@ranges@NAMES == "HXB2_VPR"]
        NT@ranges@NAMES[NT@ranges@NAMES == "NONHXB2_VPR"] <- "HXB2_VPR"
      }
      
      NT_df <- as.data.frame(NT[1])
      
      ##extract sequence from original fasta
      fasta_raw <- as.matrix(NT)
      fasta_raw <- t(fasta_raw)
      fasta_raw <- cbind(1:nrow(fasta_raw), fasta_raw, 1:nrow(fasta_raw))
      colnames(fasta_raw)[1] <- "fasta_indx"
      colnames(fasta_raw)[3] <- "freq_indx"
      
      NT_freq_cov <- cbind(rep(fasta_raw[, 2], each = 5), 
                           as.data.frame(rep(c("A", "C", "G", "T", "-"), nrow(fasta_raw)))
                           )
      colnames(NT_freq_cov) <- c("true_variant", "variant")

      NT_freq_cov <- NT_freq_cov %>%
        mutate(
          position = sort(rep(1:nrow(fasta_raw), 5)),
          depth = 10000,
          reference = "K03455.1"
        ) %>%
          mutate(frequency = ifelse(variant == true_variant, 1, 0)) %>%
          mutate(freq_indx = cur_group_id(), .by = position) %>%
          dplyr::select(position,
                      variant,
                      frequency,
                      reference,
                      depth,
                      freq_indx)
      
    }
    
    
    NT_tab <- NT_df %>%
      separate(x, 
               into = paste("V", 1:(max(nchar(NT_df$x))), sep = ""), 
               sep = "")
    
    NT_tab <- as.data.frame(t(NT_tab))
    
    NT_tab$row_num <- seq(1:nrow(NT_tab))
    
    NT_tab$index_fasta_raw <- NA
    NT_tab$index_fasta_raw[NT_tab[, 1] != "-" &
                             NT_tab[, 1] != "!"] <- fasta_raw[, 1] %>% unlist()
    NT_tab$index_fasta_freq = NA
    NT_tab$index_fasta_freq[NT_tab[, 1] != "-" &
                              NT_tab[, 1] != "!"] <- fasta_raw[, 3] %>% unlist()
    
    ##AA single seqs
    if (uuid_name != paste0("HXB2_", toupper(region))) {
      aas <- Biostrings::readBStringSet(
        paste0(
          "codon_align/codon_align_",
          region,"_AA/",
          uuid_name,
          ".fa"
        )
      )
      aas <- aas[!grepl(region, aas@ranges@NAMES, ignore.case = TRUE)]
      aas_df <- as.data.frame(aas[1])
    }
    if (uuid_name == paste0("HXB2_", toupper(region))) {
      #random sequence to get an hxb2 protein sequence, from the next one or otherwise the one before
      aas <- Biostrings::readAAStringSet(
        paste0(
          "references_seqs/hxb2_",
          region,
          "_aa.fa"
        )
      )
      
      ##for vpr due to frame shift, maybe hxb2 reference is not necessary at all, but lets remove it here and take only the non-frameshift complete vpr sequence
      if (region == "vpr") {
        aas <- aas[!aas@ranges@NAMES == "HXB2_VPR"]
        aas@ranges@NAMES[aas@ranges@NAMES == "NONHXB2_VPR"] <- "HXB2_VPR"
      }
      
      aas_df <- as.data.frame(aas[1])
    }
    
    aas_tab <- aas_df %>% 
      separate(x, into = paste("V", 1:(max(nchar(aas_df$x))), sep = ""), sep = "") 
    
    aas_tab <- as.data.frame(t(aas_tab))
    
    aas_tab$row_num <- seq(1:nrow(aas_tab))
    
    aas_tab <- aas_tab[rep(seq_len(nrow(aas_tab)), each = 3), ]
    
    aas_tab$row_num_nt <- seq(1:nrow(aas_tab))
    
    nt_aa <- cbind(aas_tab, NT_tab)
    
    colnames(nt_aa) <- c(
      "id_aa",
      "row_num_aa",
      "row_num_nt",
      "id_nt",
      "row_num_ignore",
      "index_fasta_raw",
      "index_fasta_raw_freq"
    )
    
    ##combine all
    nt_aa$row_num_aa_al <- NA
    orgindex <- 0
    for (i in 1:nrow(aa_tab)) {
      amino = aa_tab[i, 1]
      
      orgindex_temp = which(nt_aa$id_aa == amino)
      orgindex_temp = orgindex_temp[orgindex_temp > orgindex]
      orgindex = min(orgindex_temp)
      
      nt_aa$row_num_aa_al[orgindex:(orgindex + 2)] =  aa_tab[i, 2]
      
      orgindex = orgindex + 2
    }
    
    nt_aa = nt_aa %>% 
      mutate(nt_al = ifelse(is.na(row_num_aa_al), "-", id_nt))
    
    nt_aa = nt_aa %>% 
      filter(nt_al != "-") %>% 
      dplyr::select(row_num_aa_al,
      nt_al,
      index_fasta_raw,
      index_fasta_raw_freq)
    
    nt_aa$rowindx = seq(1:nrow(nt_aa)) %% 3
    nt_aa$rowindx[nt_aa$rowindx == 0] = 3
    
    nt_aa$nt_indx = nt_aa$row_num_aa_al * 3 + nt_aa$rowindx - 3
    
    ref_seq_frame = as.data.frame(seq(1:6000))
    colnames(ref_seq_frame) = "ref_frame"
    
    ref_seq_frame = merge(
      ref_seq_frame,
      nt_aa,
      by.x = "ref_frame",
      by.y = "nt_indx",
      all = TRUE) %>%
      mutate(nt_al = ifelse(is.na(nt_al), "-", nt_al)) %>%
      mutate(index_fasta_raw_freq = as.numeric(index_fasta_raw_freq))
    
    nt_al_list <- append(nt_al_list, paste0(ref_seq_frame$nt_al, collapse = ""))
    seq_ids <- c(seq_ids, uuid_name)
    ##format frequency format gwas file
    aligned_freq <- inner_join(
      NT_freq_cov,
      ref_seq_frame,
      by = c("freq_indx" = "index_fasta_raw_freq")
    )
    ##from long to wide format
    aligned_freq <- aligned_freq %>% 
      dplyr::select(ref_frame, variant, frequency, depth) %>% 
      pivot_wider(
        names_from = c(ref_frame, variant),
        values_from = c(frequency, depth)
      )
    ##extract depth columns, for reduction from 5 to one depth column
    depth_ps <- aligned_freq %>% select(contains("depth"))
    ##take only one for each position, I just choose "G". (could be any, depth is the same for all)
    depth_ps <- depth_ps %>% select(contains("G"))
    colnames(depth_ps) = gsub("_G", "", colnames(depth_ps))
    ##remove orginal depths and add formatted ones
    aligned_freq <- aligned_freq %>% select(-contains("depth"))
    print("5")
    ##combine sequence name, frequencies, and sequencing depth
    aligned_freq <- cbind(uuid_name, aligned_freq, depth_ps)
    print("6")
    ##combine all sequences
    freq_aligned_all <- bind_rows(freq_aligned_all, aligned_freq)
    
    print(uuid_name)
  }
  
  file_name_output <- gsub("alignment_aa", "alignment_basefreq", file_name)
  file_name_output <- gsub(".fa", ".csv", file_name_output)
  write.csv(freq_aligned_all, paste0(gsub("alignment_aa","alignment_basefreq",outputpath), "/", file_name_output))
  
  seqinr::write.fasta(
    as.list(nt_al_list),
    names = c(seq_ids),
    as.string = FALSE,
    file.out = paste0(gsub("alignment_aa","alignment_nt",outputpath), "/", file_name)
  )
  
  ##decipher post msa processing
  print("start Decipher alignment adjustment")
  dna <- readDNAStringSet(paste0(gsub("alignment_aa","alignment_nt",outputpath), "/",file_name))
  
  for (l in 1:10) {
    dna <- DECIPHER::AdjustAlignment(dna)
  }
  
  writeXStringSet(dna, paste0(gsub("alignment_aa","alignment_nt",outputpath), "/", file_name))
}

# ---------------------------
# 3. Run aa to nt
# ---------------------------

AAtoNT(region = region,
       depth_threshold = 20)

print(sessionInfo())