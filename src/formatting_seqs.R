############################################################
# Project: SHCGWAS
# Script:  formatting_seqs.R
# Author:  Marius Zeeb
# Date:    2025-09-22
# Purpose: Generates fasta files from frequency files from NGS assemblies
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

raw_sequences_path <- args[1]      
freq_sequences_path <- args[2]  
fasta_sequences_path <- args[3]  

print(args)

# Load libraries
library(foreach)
library(doParallel)
library(tidytable)

registerDoParallel(cores = 8L)

# ---------------------------
# 2. Generate frequency files
# ---------------------------

##List all raw files 
raw_seq_files <- list.files(raw_sequences_path, full.names = TRUE)

##warnings not relevant
foreach(freq_file = raw_seq_files) %dopar% {

  marytesat <- read.csv(freq_file)
  if(is.na(colnames(marytesat)[2]) || colnames(marytesat)[2] != "POS"){
    why <- why + 1
    return(NA)
  }

  marytesat <- marytesat %>%
    separate(ALT, c("first", "second","third","fourth"), sep = "/", remove = FALSE) %>%
    separate(AF, c("first_f", "second_f","third_f","fourth_f"), sep = "/", remove = FALSE) %>%
    mutate(first_f = as.numeric(first_f),
           second_f = as.numeric(second_f),
           third_f = as.numeric(third_f),
           fourth_f = as.numeric(fourth_f)) %>%
    mutate(ref_f = 100 - rowSums(dplyr::select(., first_f, second_f, third_f,fourth_f), na.rm = TRUE)) %>%
    mutate(A = case_when(REF == "A" ~ paste0(ref_f),
                         first == "A" ~ paste0(first_f), 
                         second == "A" ~ paste0(second_f),
                         third == "A" ~ paste0(third_f),
                         fourth == "A" ~ paste0(fourth_f))) %>%
    mutate(T = case_when(REF == "T" ~ paste0(ref_f),
                         first == "T" ~ paste0(first_f), 
                         second == "T" ~ paste0(second_f),
                         third == "T" ~ paste0(third_f),
                         fourth == "T" ~ paste0(fourth_f))) %>%
    mutate(G = case_when(REF == "G" ~ paste0(ref_f),
                         first == "G" ~ paste0(first_f), 
                         second == "G" ~ paste0(second_f),
                         third == "G" ~ paste0(third_f),
                         fourth == "G" ~ paste0(fourth_f))) %>%
    mutate(C = case_when(REF == "C" ~ paste0(ref_f),
                         first == "C" ~ paste0(first_f), 
                         second == "C" ~ paste0(second_f),
                         third == "C" ~ paste0(third_f),
                         fourth == "C" ~ paste0(fourth_f))) %>%
    mutate(GAP = NA) %>%
    dplyr::select(POS,A,C,G,T,GAP,COV)

  marytesat <- marytesat %>%
    pivot_longer(., cols = c("A", "C", "G", "T", "GAP") , names_to = "variable", values_to = "value") %>%
    arrange(POS, match(variable, c("A", "C", "G","T", "GAP")))

  marytesat$variable <- as.character(marytesat$variable)
  marytesat$variable[which(marytesat$variable == "GAP")] <- "-"
  marytesat$value[which(is.na(marytesat$value))] <- 0
  marytesat$value <- as.numeric(marytesat$value) / 100
  marytesat$COV[which(is.na(marytesat$COV))] <- 0
  marytesat$reference <- "K03455.1"

  freq <- marytesat %>%
    dplyr::select(POS,variable,value,reference,COV)

  colnames(freq) <- c("position", "variant", "frequency", "reference", "depth")

  output_path <- gsub(raw_sequences_path, freq_sequences_path, freq_file)
  write.csv(freq, output_path, row.names = FALSE)

}

# ---------------------------
# 3. Generate fasta files
# ---------------------------

##List all frequency files 
freq_seq_files <- list.files(freq_sequences_path, full.names = TRUE)

foreach(fasta_file = freq_seq_files) %dopar% {

  ##load frequency file
  cte <- read.csv(file = fasta_file)
  cte <- cte %>%
    arrange(position, desc(frequency)) %>%
    slice(1, .by = position) %>%
    mutate(variant = replace(variant, depth <= 20, "N")) %>% #coverage filter (alternatively just remove)
    dplyr::select(variant)

  fasta_internal_name <- gsub(paste0(freq_sequences_path,"|/|.csv"),"", fasta_file)

  fasta <- paste0(">",fasta_internal_name,"\n",
                 paste(cte$variant, collapse = ""))

  output_path <- gsub(freq_sequences_path, fasta_sequences_path, fasta_file)
  output_path <- gsub(".csv", ".fa", output_path)
  write.table(fasta,
              output_path,
              quote = FALSE,
              row.names = FALSE,
              col.names  = FALSE)
}

print(sessionInfo())

############################################################
# End of Script
############################################################