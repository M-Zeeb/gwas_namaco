############################################################
# Project: SHCGWAS
# Script:  apd_score.R
# Author:  Marius Zeeb
# Date:    2025-10-31
# Purpose: Calculates APD score
############################################################

# ---------------------------
# 1. Setup
# ---------------------------
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

freq_sequences_path <- args[1]
fasta_sequences_path <- args[2]

# Load libraries
library(devtools)
library(foreach)
library(doParallel)
library(Biostrings, include.only = c("xscat", "writeXStringSet", "readDNAStringSet"))
library(tidytable)

registerDoParallel(cores = 7L)

# ---------------------------
# 2. Functions
# ---------------------------

## APD tiydtable version ----
apd_calc = function(uuid,region,depth_threshold){
  #uuid = "0031606e-edb4-45c0-94cc-a7b793153c70"
  #region = "pol"
  #depth_threshold = 100

  ##count third codon position
  co_al = as.data.frame(
    read.table(
      paste0("codon_align/",
             "codon_align_",region,"_NT/",uuid,".fa"))[4,])
  colnames(co_al) = "seq"
  co_al = co_al %>% separate(seq, into = paste("POS", 0:(max(nchar(co_al$seq))), sep = ""), sep = "") %>% 
    #select(-'POS0') %>% 
    tidyr::gather(number, base) %>% 
    mutate(number = rownames(.))
  co_al <- co_al %>% 
    filter(base %in% c("A","G","C","T","!","N")) %>%
    ##probably would be more elegant with modulo operator or something
    mutate(co_pos = rep(seq(1:3),nrow(.)/3)) %>%
    mutate(pos_codon = paste0(number,"_",co_pos)) %>%
    filter(base %in% c("A","G","C","T","N"))

  ##get originial assembly position
  seq_org_pos <- readDNAStringSet(paste0(fasta_sequences_path,uuid,".fa"))
  seq_meta <- read.csv(paste0("blastmeta/blastmeta_",region,".csv")) %>% 
    filter(QueryID == uuid)
  seq_org_pos <- as.data.frame(chartr("N", "N", seq_org_pos))
  colnames(seq_org_pos) <- "seq"
  seq_org_pos <- seq_org_pos %>% 
    separate(seq, into = paste("POS", 0:(max(nchar(seq_org_pos$seq))), sep = ""), sep = "") %>% 
    #select(-'POS0') %>% 
    tidyr::gather(number_assembly, base_assembly) %>% 
    mutate(number_assembly =  rownames(.)) %>%
    filter(number_assembly %in% seq_meta$Q.start:seq_meta$Q.end)

  ambg_pos <- as.data.frame(cbind(co_al,seq_org_pos)) %>% 
    filter(co_pos == 3) %>% 
    filter(base != "N")

  ##load assembly information
  freqs <- read.csv(paste0(freq_sequences_path,uuid,".csv")) %>%
    filter(depth >= depth_threshold)

  ##calculate score
  ambg_score <- freqs %>% 
    filter(position %in% ambg_pos$number_assembly) %>%
    mutate(countit = max(frequency) < 0.99, .by = position) %>%
    mutate(score_base = frequency * (1-frequency)) %>%
    mutate(ambg_score_pos = sum(score_base), .by = position) %>%
    mutate(ambg_score_pos = ifelse(countit == FALSE,0,ambg_score_pos)) %>%
    tidytable::slice(1, .by = position) %>%
    mutate(n = n()) %>%
    mutate(apd_score = sum(ambg_score_pos)/n()) %>%
    tidytable::slice(1)

  return(c(ambg_score$apd_score,ambg_score$n))
}

# ---------------------------
# 3. Calculate APDs for pol, env, and gag
# ---------------------------

uuids_pol <- read.csv(paste0("blastmeta/blastmeta_pol.csv")) %>% 
  mutate(region = "pol") %>% filter(!is.na(SubjectID))
uuids_gag <- read.csv(paste0("blastmeta/blastmeta_gag.csv")) %>%
  mutate(region = "gag") %>% filter(!is.na(SubjectID))
uuids_env <- read.csv(paste0("blastmeta/blastmeta_env.csv")) %>%
  mutate(region = "env") %>% filter(!is.na(SubjectID))

uuid_apd <- rbind(
      uuids_pol,
      uuids_gag,
      uuids_env)

apds <- foreach(x = 1:nrow(uuid_apd), .combine = "rbind") %dopar% {

  apd <- apd_calc(uuid = uuid_apd$QueryID[x],region = uuid_apd$region[x], depth_threshold = 100)
  if(length(apd) == 0){ apd <- c("fail","fail")}

  return(cbind(uuid_apd$QueryID[x],uuid_apd$region[x],apd))

}

apds_format <- apds %>% 
  as.data.frame() %>% 
  mutate(n = lead(apd))
apds_format <- apds_format[seq(1,nrow(apds_format), by = 2),]
apds_format <- apds_format %>% 
  pivot_wider(names_from = V2, values_from = c(apd,n))

colnames(uuids_pol) = paste0("pol_",colnames(uuids_pol))
colnames(uuids_gag) = paste0("gag_",colnames(uuids_gag))
colnames(uuids_env) = paste0("env_",colnames(uuids_env))

apds_format = apds_format %>%
  merge(.,uuids_pol, by.x = "V1",by.y = "pol_QueryID", all.x = TRUE, all.y = FALSE) %>%
  merge(.,uuids_gag, by.x = "V1",by.y = "gag_QueryID", all.x = TRUE, all.y = FALSE) %>%
  merge(.,uuids_env, by.x = "V1",by.y = "env_QueryID", all.x = TRUE, all.y = FALSE)
colnames(apds_format)[1] = "uuid"

write.csv(apds_format,
          paste0("results/apd_scores.csv")
          )

print(sessionInfo())