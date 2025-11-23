############################################################
# Project: SHCGWAS
# Script:  plink.R
# Author:  Marius Zeeb
# Date:    2025-09-25
# Purpose: Performs sequence formatting with plink for PCA calculation
############################################################

# ---------------------------
# 1. Setup
# ---------------------------

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

region <- args[1]  
res <- NULL

library(tidyverse)
library(svMisc)
library(foreach)
library(doParallel)
library(Rcpp)
library(qdapTools)
library(glmnet)
library(fastDummies)

#sourceCpp("/Users/mariuszeeb/Downloads/bacterial-heritability-main/tableC.cpp")

# ---------------------------
# 2. Functions
# ---------------------------

##reformat sequence datas

reformat_seq_to_csv <- function(msa){
  
  msa_df <- as.data.frame(msa)
  msa_df$V1 <- msa@ranges@NAMES
  ##Rename column names
  names(msa_df) <- c("seq","id")
  msa_df <- msa_df %>% 
    dplyr::select(id,seq)
  ##Dataset for sequencing data
  msa <- msa_df %>% 
    separate(seq, into = paste("V", 1:(max(nchar(msa_df$seq))), sep = ""), sep = "")
  
  return(msa)
}

read_alignment_and_format <- function(region){
  
  ##load alignment
  msa <- Biostrings::readAAStringSet(paste0("alignment_nt/alignment_",region,".fa"))
  
  ##reformat dnastring-format to dataframe
  msa <- reformat_seq_to_csv(msa)
  write.csv(msa,paste0("alignment_nt/alignment_",region,".csv"))
  
  return(msa)
}

##store baseposition in relation to reference for later seqs to 
ref_seq_pos <- function(msa,region){
  
  hxb2_region <- paste0("HXB2_",toupper(region))
  
  seqs_to_keep_neuro <- as.matrix(msa %>%
                                   filter(id == hxb2_region))
  
  seqs_to_keep_neuro <- t(seqs_to_keep_neuro)[-1,,drop = FALSE] %>% as.data.frame()
  colnames(seqs_to_keep_neuro) <- "variant"
  seqs_to_keep_neuro$gwas_pos <- as.numeric(gsub("V", "", rownames(seqs_to_keep_neuro)))  
  
  seqs_to_keep_neuro  <- seqs_to_keep_neuro  %>%
    filter(variant != "-") %>%
    mutate(hxb2_index = 1:n()) %>%
    mutate(hxb2_index_overall = as.numeric(hxb2_index))
  
  seqs_to_keep_neuro$msa_pos <- seqs_to_keep_neuro$gwas_pos 
  write.csv(seqs_to_keep_neuro,paste0("alignment_nt/seqs_to_keep_",region,".csv"))
  
  return(seqs_to_keep_neuro)
}

##resitances definitions for pol only
resistances_pos <- function(region){
  
  ##resistance position
  resistances <- as.data.frame(c(89,90,91,95,96,97,98,99,100,137,138,139,140,141,142,143,144,145,149,
                                150,151,161,162,163,173,174,175,221,222,223,227,228,229,245,246,247,
                                248,249,250,251,252,253,263,264,265,269,270,271,419,420,421,482,483,
                                484,491,492,493,497,498,499,506,507,508,518,519,520,521,522,523,527,
                                528,529,596,597,598,599,600,601,605,606,607,614,615,616,620,621,622,
                                641,642,643,644,645,646,710,711,712,749,750,751,833,834,835,839,840,
                                841,848,849,850,860,861,862,866,867,868,926,927,928,941,942,943,953,
                                954,955,959,960,961,971,972,973,977,978,979,986,987,988))
  colnames(resistances) <- "respos"
  if(region == "pol"){
    resistances = resistances %>% 
      mutate(respos3 = respos * 3) %>% 
      mutate(respos2 = respos3 - 1) %>% 
      mutate(respos1 = respos3 - 2)
  }

  if(region == "whole"){
    resistances <- resistances %>% 
      mutate(respos3 = (respos * 3) + 2083) %>% 
      mutate(respos2 = respos3 - 1) %>% 
      mutate(respos1 = respos3 - 2)
  }
  
  return(resistances)
}

##plinkformatting
plink_formatting_vcf_alleles <- function(msa,n_snps,n_samp){
  
  msa <- as.data.frame(msa) #enforce a data frame
  
  vcf_allels <- NULL
  vcf_allels <- as.data.frame(matrix(ncol = 1, nrow = n_snps-1))
  for(i in 2:n_snps){
    vcf_allels$V1[i-1] <- i-1
    vcf_allels$A[i-1] <- sum(msa[,i] == "A",na.rm = TRUE)
    vcf_allels$C[i-1] <- sum(msa[,i] == "C",na.rm = TRUE)
    vcf_allels$G[i-1] <- sum(msa[,i] == "G",na.rm = TRUE)
    vcf_allels$T[i-1] <- sum(msa[,i] == "T",na.rm = TRUE)
    vcf_allels$N[i-1] <- sum(is.na(msa[,i]) | (!msa[,i] %in% c("A","G","C","T")))
  }
  
  return(vcf_allels)
}

generate_vcf_empty <- function(vcf_allels,maf,n_snps,n_samp){
  
  ##vcf file
  gwas_vcf <- NULL
  gwas_vcf <- as.data.frame(matrix(ncol = (n_samp+9), nrow = (n_snps-1)))
  for(z in 1:(n_snps-1)){
    gwas_vcf[z,1] <- 1
    gwas_vcf[z,2] <- z 
    gwas_vcf[z,3] <- "." 
    gwas_vcf[z,4] <- colnames(vcf_allels)[which(vcf_allels[z,2:5] == max(vcf_allels[z,2:5]))[1]+1]
    gwas_vcf[z,5] <- paste0(unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][1]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][2]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][3]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][4]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][5]))
    gwas_vcf[z,5] <- gsub(pattern = paste0("NA,|NA|,NA"),"",gwas_vcf[z,5])
    gwas_vcf[z,5] <- gsub(pattern = paste0("\\",gwas_vcf[z,4],","),"",gwas_vcf[z,5])
    gwas_vcf[z,5] <- gsub(pattern = paste0(",\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
    gwas_vcf[z,5] <- gsub(pattern = paste0("\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
    gwas_vcf[z,5][gwas_vcf[z,5] == ""] <- "."
    gwas_vcf[z,6] <- "."
    gwas_vcf[z,7] <- "."
    gwas_vcf[z,8] <- "."
    gwas_vcf[z,9] <- "GT"
    
  }  
  colnames(gwas_vcf)[1:9] = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  
  return(gwas_vcf)
}

generate_vcf_fill <- function(msa,gwas_vcf,n_snps,n_samp){
  
  for(z in 1:nrow(msa)){
    colnames(gwas_vcf)[z+9] <- as.character(msa$id[z])
    chrs <- strsplit(as.character(paste0(gwas_vcf[1:(n_snps-1),4],",",gwas_vcf[1:(n_snps-1),5],",",msa[z,2:n_snps])), split = ",")
    gwas_vcf[,z+9] <- paste0(unlist(lapply(chrs,anyDuplicated, fromLast = TRUE))-1)
    
    progress(z,progress.bar = FALSE, init = (z == 1), max.value = n_samp)
  }
  
  gwas_vcf[,][gwas_vcf[,] == "-1"] <- "." 
  
  return(gwas_vcf) 
}

binarising <- function(gwas_vcf,altc,n_snps,n_samp){
  
  othercs <- c(1,2,3,4)[c(1,2,3,4) != altc] 
  gwas_vcf_alt2 <- gwas_vcf
  nseq <- n_samp + 9 
  gwas_vcf_alt2$ID <- NA
  gwas_vcf_alt2$ID[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)] <- paste0(gwas_vcf_alt2$POS[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],
                                                                                                        gwas_vcf_alt2$REF[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],
                                                                                                        altc,
                                                                                                        sapply(
                                                                                                          strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),
                                                                                                          "[[", 
                                                                                                          altc)
                                                                                                        ) 
  
  gwas_vcf_alt2 <- gwas_vcf_alt2[which(!is.na(gwas_vcf_alt2$ID)),]
  
  ##to keep non-missing with other base as 0 and complete missing as .
  truemiss <- stringr::str_locate(pattern = "N", gsub(",","",gwas_vcf_alt2[,"ALT"]))[,1]
  truemiss[is.na(truemiss)] <- 99
  miss <- gwas_vcf_alt2[,10:nseq] == truemiss
  
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[1] | miss)] <- "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[1] | (!miss))] <- "0" 
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[2] | miss)] <- "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[2] | (!miss))] <- "0"
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[3] | miss)] <- "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[3] | (!miss))] <- "0"
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == altc] <- 1
  gwas_vcf_alt2$ALT <- sapply(strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),"[[", altc)
  
  return(gwas_vcf_alt2)
  
}

final_vcf <- function(gwas_vcf_alt,seqs_to_keep_neuro_todel,region,n_samp,res=NULL){
  
  gwas_vcf_alt$`#CHROM` <- 1
  ##remove resistance position
  if((!is.null(res)) & ((region == "pol") | (region == "whole"))){
    gwas_vcf_alt <- gwas_vcf_alt %>% filter(POS %in% seqs_to_keep_neuro_todel$msa_pos)
  }
  gwas_vcf_alt <- gwas_vcf_alt %>% 
    filter(ALT != "N" & ALT != ".")
  gwas_vcf_alt$FILTER <- "PASS"
  gwas_vcf_alt$POS <- c(1:nrow(gwas_vcf_alt))
  
  vcf_body_path <- paste0("pca/pca_",region,"/",
                    paste(c("plink2vcf",region,res), collapse = "_"),
                    ".vcf")
  
  write.table(gwas_vcf_alt,vcf_body_path,
              sep ="\t", quote = FALSE, row.names = FALSE)
  
  header <- paste0("##fileformat=VCFv4.2 
##fileDate=",Sys.Date(),"
##source=",version$version.string,"
##contig=<ID=1,length=",n_samp,">
##INFO=<ID=PR,Number=0,Type=Flag,Description=\"",paste(c(region,res),collapse = "_"),"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
  
  vcf_header_path <- paste0("pca/pca_",region,"/",
                           paste(c("vcfheader",region,res), collapse = "_"),
                           ".txt")
  
  write.table(header,vcf_header_path,
              sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  vcf_path <- paste0("pca/pca_",region,"/",
                           paste(c("vcf",region,res), collapse = "_"),
                           ".vcf")
  
  ##merge vcf with fileheader
  system(paste0("cat ",vcf_header_path," ",vcf_body_path," > ",vcf_path))
  
}

final_vcf_x <- function(gwas_vcf_alt,seqs_to_keep_neuro_todel,region,n_samp,res=NULL){

  gwas_vcf_alt$`#CHROM` <- 1
  ##remove resistance position
  if((!is.null(res)) & ((region == "pol") | (region == "whole"))){
    gwas_vcf_alt <- gwas_vcf_alt %>% filter(POS %in% seqs_to_keep_neuro_todel$msa_pos)
  }
  gwas_vcf_alt <- gwas_vcf_alt %>% 
    filter(ALT != "N" & ALT != ".")
  gwas_vcf_alt$FILTER <- "PASS"
  gwas_vcf_alt$POS <- c(1:nrow(gwas_vcf_alt))
  gwas_vcf_alt[,10:ncol(gwas_vcf_alt)][gwas_vcf_alt[,10:ncol(gwas_vcf_alt)] == 1] <- 2
  
  vcf_body_path <- paste0("pca/pca_",region,"/",
                         paste(c("plink2vcf",region,res), collapse = "_"),
                         "_x.vcf")
  write.table(gwas_vcf_alt,vcf_body_path,
              sep ="\t", quote = FALSE, row.names = FALSE)
  
  header <- paste0("##fileformat=VCFv4.2 
##fileDate=",Sys.Date(),"
##source=",version$version.string,"
##contig=<ID=1,length=",n_samp,">
##INFO=<ID=PR,Number=0,Type=Flag,Description=\"",paste(c(region,res),collapse = "_"),"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
  
  vcf_header_path <- paste0("pca/pca_",region,"/",
                           paste(c("vcfheader",region,res), collapse = "_"),
                           "_x.txt")
  
  write.table(header,vcf_header_path,
              sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ##merge vcf with fileheader
  vcf_path <- paste0("pca/pca_",region,"/",
                    paste(c("vcf",region,res), collapse = "_"),
                    "_x.vcf")
  system(paste0("cat ",vcf_header_path," ",vcf_body_path," > ",vcf_path))
  
}

sex_update <- function(gwas_vcf_alt,region){
  
  sex_data <- as.data.frame(colnames(gwas_vcf_alt)[10:ncol(gwas_vcf_alt)])
  sex_data <- sex_data %>% 
    mutate(V1 = 0, V3 = 1) %>% 
    dplyr::select(c(2,1,3))
  
  write.table(sex_data,
              paste0("pca/pca_",region,"/",
                     paste(c("sex_update","gcta",region,res),collapse = "_")),
              sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

plink_dummy_phenotype <- function(gwas_vcf_alt,region){
  
  dummy_data <- as.data.frame(colnames(gwas_vcf_alt)[11:ncol(gwas_vcf_alt)])
  dummy_data <- dummy_data %>% 
    mutate(dummy_pheno = 3) %>% 
    dplyr::select(c(1,2))
  colnames(dummy_data) <- c("#IID","dummy_pheno")
  
  write.table(dummy_data,
              paste0("pca/pca_",region,"/",
                     paste(c("dummy_pheno",region),collapse = "_"),
                     ".txt"),
              sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}

##region possibilities ("tat","pol","env","gag","whole"), res "_res" for pol, type ("_whole" or "_resdb")
main <- function(region,res = NULL){
  
  print(paste0(region," start"))
  
  ##read and format
  msa <- read_alignment_and_format(region)

  ##store ref seq pos
  seqs_to_keep <- ref_seq_pos(msa,region)
  
  ##when region is pol
  seqs_to_keep_neuro_todel <- NULL
  if((region == "pol") | (region == "whole")){
    
    ##resistance definitions
    resistances <- resistances_pos(region)
    
    seqs_to_keep_neuro_todel <- seqs_to_keep %>% 
      filter(!(hxb2_index %in% resistances$respos1 | 
                 hxb2_index %in% resistances$respos2 | 
                 hxb2_index %in% resistances$respos3))
  }
  
  ##number of snps
  n_snps <- length(msa)
  #number of samplesß
  n_samp <- nrow(msa)
  
  ##ambguious nucleotides as "N"
  msa <- msa %>% 
    mutate(across(2:(n_snps), ~ replace(., which(is.na(.) | !. %in% c("A","C","G","T")),"N")))
  
  ##make table with nucleotides at each position across sequence samples 
  vcf_alleles <- plink_formatting_vcf_alleles(msa,n_snps = n_snps,n_samp = n_samp)
  
  ##create empty vcf file
  ##maf is minor allele frequency threshold
  gwas_vcf <- generate_vcf_empty(vcf_alleles,maf=5,n_snps = n_snps,n_samp = n_samp)
  
  ##file empty vcf file
  gwas_vcf <- generate_vcf_fill(msa,gwas_vcf,n_snps = n_snps,n_samp = n_samp)
  
  ##binarise vcf file
  ##altc for different possible minor alleles
  gwas_vcf_bin = rbind(binarising(gwas_vcf = gwas_vcf, altc = 1,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 2,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 3,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 4,n_snps = n_snps,n_samp = n_samp))   
  
  ##generate final vcf file for processing with plink
  final_vcf(gwas_vcf_bin,seqs_to_keep_neuro_todel,region,n_samp)
  ##for weird plink behaviour
  final_vcf_x(gwas_vcf_bin,seqs_to_keep_neuro_todel,region,n_samp)
  
  ##sex update
  sex_update(gwas_vcf_bin,region)
  
  ##dummy phenotype
  plink_dummy_phenotype(gwas_vcf_bin,region)

  
  print(paste0(region," plink preformatting done"))
  
}

# ---------------------------
# 3. Run main
# ---------------------------

main(region = region, res = res)

print(sessionInfo())