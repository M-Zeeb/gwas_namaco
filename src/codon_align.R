############################################################
# Project: SHCGWAS
# Script:  codon_align.R
# Author:  Marius Zeeb
# Date:    2025-09-22
# Purpose: Generates AA alignments from fasta files
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

registerDoParallel(cores = 8L)

# ---------------------------
# 2. Functions
# ---------------------------

##blast----
region_extract = function(uuids,region,ref_data_base) {
   #region = "pol"
   #ref_data_base = "hxx"
   #uuids = uuids$uuid
  
  outputpath <- paste0("blasted_seqs/blasted_seqs_",region,"/")
  
  ##BLAST database
  rBLAST::makeblastdb(paste0("references_seqs/",ref_data_base,"_",region,".fa"),dbtype = "nucl")
  bl <- rBLAST::blast(paste0("references_seqs/",ref_data_base,"_",region,".fa"), type = "blastn")
  blastmeta = data.frame(QueryID = character(),
                         SubjectID = character(),  
                         Perc.Ident = numeric(),  
                               Alignment.Length = integer(),
                               Mismatches = integer(),   
                               Gap.Openings = integer(),  
                               Q.start  = integer(),      
                               Q.end   = integer(),        
                               S.start  = integer(),        
                               S.end   = integer(),         
                               E = numeric(),             
                               Bits = numeric(), 
                               stringsAsFactors = FALSE)
  y = 1
  blastmeta = foreach(i = uuids, .combine = "rbind") %dopar%{

    ##sequences
    seq <- Biostrings::readDNAStringSet(paste0(fasta_sequences_path,i,".fa"))
    seq_df = as.data.frame(Biostrings::readDNAStringSet(paste0(fasta_sequences_path,i,".fa")))

    ##BLAST

    blast_seq = tryCatch({predict(bl, seq, BLAST_args = "-max_target_seqs 2000 -evalue 4000")}, 
             error = function(e) {
                print("why")
                })

    if("X1" %in% colnames(blast_seq)){
      colnames(blast_seq) <- c("QueryID","SubjectID","Perc.Ident","Alignment.Length","Mismatches","Gap.Openings","Q.start","Q.end","S.start","S.end","E","Bits")
    }

    if(!"Alignment.Length" %in% colnames(blast_seq)){
      blast_seq <- blast_seq %>%
        rename(Alignment.Length = length,
               E = evalue,
               QueryID = qseqid,
               SubjectID = sseqid,
               Perc.Ident = pident,
               Mismatches = mismatch,
               Gap.Openings = gapopen,
               Q.start = qstart,
               Q.end = qend,
               S.start = sstart,
               S.end = send,
               Bits = bitscore)
    }

    if(nrow(blast_seq) == 0){
      blast_seq[nrow(blast_seq) + 1, ] <- NA
      blast_seq[1,1] = i
      #next
      return(blast_seq)
    }

    ##seems to be the most practical to sort for length and then for E score
    blast_seq <- blast_seq %>% arrange(desc(Alignment.Length),E) %>% slice(1)
    #blast_seq$Q.end[1] = 5101
    ##MERGE BLAST with raw sequence
    seq_df %>%
      mutate(reg_seq = substr(x, blast_seq$Q.start, blast_seq$Q.end))  %>%
      mutate(forcsv = paste0(">",rownames(seq_df)[1],"\n",reg_seq)) %>% {
      write.table(.$forcsv,paste0(outputpath,i,".fa"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
      }
    #blastmeta = rbind(blastmeta,blast_seq)
    return(blast_seq)
  }
  write.csv(blastmeta,
            paste0(paste(c("blastmeta/blastmeta",region),collapse = "_"),".csv"),
            row.names = FALSE)
}

##exon merge----
exon_merge = function(uuids,region){

  blastmeta_exon1_path <- paste(c("blastmeta/blastmeta",tolower(region),"exon1"),collapse = "_")
  blastmeta_exon2_path <- paste(c("blastmeta/blastmeta",tolower(region),"exon2"),collapse = "_")

  blastseq_exon1_path <- paste(c("blasted_seqs/blasted_seqs",tolower(region),"exon1"), collapse = "_")
  blastseq_exon2_path <- paste(c("blasted_seqs/blasted_seqs",tolower(region),"exon2"), collapse = "_")

  outputpath <- paste(c("blasted_seqs/blasted_seqs",tolower(region)), collapse = "_")

  blast_ex1 <- read.csv(paste0(blastmeta_exon1_path,".csv")) %>% filter(!is.na(SubjectID))
  blast_ex2 <- read.csv(paste0(blastmeta_exon2_path,".csv")) %>% filter(!is.na(SubjectID))

  uuid_temp <- uuids %>% as.data.frame() %>% tidytable::rename("uuid" = 1) %>%
    mutate(ex1 = uuid %in% blast_ex1$QueryID) %>%
    mutate(ex2 = uuid %in% blast_ex2$QueryID)

  for(uuid in uuid_temp$uuid[uuid_temp$ex1 == TRUE & uuid_temp$ex2 == TRUE]){
    prot_comb <- xscat(
                  readDNAStringSet(
                      paste0(blastseq_exon1_path,"/",uuid,".fa")),
                  readDNAStringSet(
                      paste0(blastseq_exon2_path,"/",uuid,".fa")))
    prot_comb@ranges@NAMES <- uuid
    writeXStringSet(prot_comb,paste0(outputpath,"/",uuid,".fa"))                     
   }

  for(uuid in uuid_temp$uuid[uuid_temp$ex1 == TRUE & uuid_temp$ex2 == FALSE]){
    prot_comb <- xscat(
                  readDNAStringSet(
                    paste0(blastseq_exon1_path,"/",uuid,".fa")))
    prot_comb@ranges@NAMES <- uuid
    writeXStringSet(prot_comb,paste0(outputpath,"/",uuid,".fa"))                        
  }

  for(uuid in uuid_temp$uuid[uuid_temp$ex1 == FALSE & uuid_temp$ex2 == TRUE]){
    prot_comb <- xscat(
                  readDNAStringSet(
                    paste0(blastseq_exon2_path,"/",uuid,".fa")))
    prot_comb@ranges@NAMES <- uuid
    writeXStringSet(prot_comb,paste0(outputpath,"/",uuid,".fa"))
  }
}

##MACSE----
codon_align = function(uuids,region){

  inputpath = paste(c("blasted_seqs/blasted_seqs",region), collapse = "_")

  outputpath = paste0("codon_align/",paste(c("codon_align",region), collapse = "_"))

  for(uuid in uuids){
    system(paste0("macse -prog alignSequences -seq references_seqs/hxb2_",region,".fa -seq_lr ",
              inputpath,"/",uuid,".fa -fs_lr 10 -stop_lr 15 -fs 100 -fs_term 100 ",
              " -out_AA ",outputpath,"_AA/",uuid,".fa",
              " -out_NT ",outputpath,"_NT/",uuid,".fa"))
  }
}

# ---------------------------
# 3. Extract sequences by regions using Blast
# ---------------------------

##define list of uuids
uuids = list.files(path = fasta_sequences_path)
uuids = gsub("\\.fa","",uuids)

##Single (internally parallelized)
##extract regions via blast (reference subtype panel) ----
region_extract(uuids,"pol","hxx")
region_extract(uuids,"env","hxx")
region_extract(uuids,"gag","hxx")
region_extract(uuids,"vif","hxx")
region_extract(uuids,"nef","hxx")
region_extract(uuids,"vpu","hxx")
region_extract(uuids,"vpr","hxx")
region_extract(uuids,"tat_exon1","hxx")
region_extract(uuids,"tat_exon2","hxx")
region_extract(uuids,"rev_exon1","hxx")
region_extract(uuids,"rev_exon2","hxx")

##merge exons from rev and tat (special case due to two codons)
exon_merge(uuids, "rev")
exon_merge(uuids, "tat")

# ---------------------------
# 4. Codon alignment with MACSE
# ---------------------------

for(region in c("env","pol","gag","vif","nef","vpu","vpr")){

  uuids = read.csv(paste0("blastmeta/blastmeta_",region,".csv")) %>%
    filter(!is.na(Alignment.Length)) %>%
    select(QueryID)

  colnames(uuids) = "uuid"

  registerDoParallel(cores = 8L)

  foreach(x = unique(uuids$uuid)) %dopar% {

    codon_align(x,region)

  }
}

##CODON alignment for rev/tat exons (special case due to two codons)
for(region in c("tat", "rev")){

    uuids1 = read.csv(paste0("blastmeta/blastmeta_",region,"_exon1",".csv")) %>% 
      filter(!is.na(Alignment.Length)) %>% select(QueryID)
    colnames(uuids1) = "uuid"
    uuids2 = read.csv(paste0("blastmeta/blastmeta_",region,"_exon2",".csv")) %>% 
      filter(!is.na(Alignment.Length)) %>% select(QueryID)
    colnames(uuids2) = "uuid"

    uuids = rbind(uuids1,uuids2)
    colnames(uuids) = "uuid"

  foreach(x = unique(uuids$uuid)) %dopar% {

    codon_align(x, region)

  }
}

print(sessionInfo())