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

## blast----
region_extract <- function(uuids, region, ref_data_base) {
   #region = "pol"
   #ref_data_base = "hxx"
   #uuids = "0b513236-5f01-49d8-9e5d-f7e8cf17f4dc"
   #fasta_sequences_path = "fasta_files/"

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
    #blast_seq %>% arrange(desc(length))

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
    blast_seq <- blast_seq %>% arrange(desc(Alignment.Length), E) %>% slice(1)
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

## minimap ----
region_extract_minimap <- function(uuids, region, ref_data_base) {
  #  region = "pol"
  #  ref_data_base = "hxx"
  #  uuids = "0fc3e13d-29e3-4b37-a487-0a924890d410"
  #  fasta_sequences_path = "fasta_files/"

  outputpath <- paste0("blasted_seqs/blasted_seqs_", region, "/")
  
  ##MMI database
  system(
    paste0(
      "minimap2 -k 8 -w 4 ",
      "-d references_seqs/", ref_data_base, "_", region, ".mmi ",
      "references_seqs/", ref_data_base, "_", region, ".fa"
    )
  )

  minimapped_meta <- foreach(i = uuids, .combine = "rbind") %dopar%{
    
    #i = "0d5a6c44-b079-480a-96a1-ee7bdc6d2b7c"
    # sequences (for later)
    seq <- Biostrings::readDNAStringSet(paste0(fasta_sequences_path, i, ".fa"))
    seq_df <- as.data.frame(Biostrings::readDNAStringSet(paste0(fasta_sequences_path, i, ".fa")))

    # minimap
    system(
      paste0(
        "minimap2 ",
        "references_seqs/", ref_data_base, "_", region, ".mmi ",
        fasta_sequences_path, i, ".fa", " > ",
        "minimapped_paf/minimapped_paf_", region, "/", i, ".paf"
      )
    )

    # read paf file
      # PAF FORMAT
      # 1 	string 	Query sequence name
      # 2 	int 	Query sequence length
      # 3 	int 	Query start (0-based; BED-like; closed)
      # 4 	int 	Query end (0-based; BED-like; open)
      # 5 	char 	Relative strand: "+" or "-"
      # 6 	string 	Target sequence name
      # 7 	int 	Target sequence length
      # 8 	int 	Target start on original strand (0-based)
      # 9 	int 	Target end on original strand (0-based)
      # 10 	int 	Number of residue matches
      # 11 	int 	Alignment block length
      # 12 	int 	Mapping quality (0-255; 255 for missing)
    tryCatch({
      mapped_seq <- read.csv(paste0("minimapped_paf/minimapped_paf_", region, "/", i, ".paf"), sep = "\t", header = FALSE)
      }, 
             error = function(e) {
                print(e$message)
                })
    
    # if no file exists or empty, create empty dataframe
    if (!exists("mapped_seq")) {
      mapped_seq <- data.frame(matrix(NA, nrow = 0, ncol = 11))
    }

    colnames(mapped_seq)[1:11] <- c("QueryID", "Q.length", "Q.start", "Q.end", "Strand", "SubjectID", "S.length", "S.start", "S.end", "ResMatches", "Alignment.Length")

    # to be compatible with blast output
    mapped_seq <- mapped_seq %>%
      mutate(
        Q.start = ifelse(Q.start > 2, Q.start - 2, Q.start - Q.start),
        Q.end = ifelse((Q.length - Q.end) > 2, Q.end + 2, Q.length),
        Alignment.Length = Q.end - Q.start,
        Perc.Ident = ResMatches / Alignment.Length * 100,
        Mismatches = Alignment.Length - ResMatches,
        Gap.Openings = NA,
        E = NA,
        Bits = NA
        ) %>%
      select(QueryID, SubjectID, Perc.Ident, Alignment.Length, Mismatches, Gap.Openings, Q.start, Q.end, S.start, S.end, E, Bits) %>%
      as.data.frame()

    mapped_seq <- mapped_seq %>%
      filter(Perc.Ident >= 15) %>%
      as.data.frame()

    if (nrow(mapped_seq) == 0) {
      mapped_seq[nrow(mapped_seq) + 1, ] <- NA
      mapped_seq[1, 1] <- i
      #next
      return(mapped_seq)
    }

    ## take mapping with longest alignment
    mapped_seq <- mapped_seq %>% 
      arrange(desc(Alignment.Length)) %>% 
      slice(1)

    ##MERGE MINIMAPPED with raw sequence
    seq_df %>%
      mutate(reg_seq = substr(x, mapped_seq$Q.start, mapped_seq$Q.end)) %>%
      mutate(forcsv = paste0(">", rownames(seq_df)[1], "\n", reg_seq)) %>% {
      write.table(.$forcsv, paste0(outputpath, i, ".fa"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
      }

    return(mapped_seq)
  }
  write.csv(minimapped_meta,
            paste0(paste(c("blastmeta/blastmeta", region), collapse = "_"),".csv"),
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
codon_align = function(uuids, region){

  inputpath = paste(c("blasted_seqs/blasted_seqs", region), collapse = "_")

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
uuids <- list.files(path = fasta_sequences_path)
uuids <- gsub("\\.fa", "", uuids)

##Single (internally parallelized)
##extract regions via blast (reference subtype panel) ----
# region_extract_minimap() with, more robust
# region_extract() with blastn
region_extract_minimap(uuids, "pol", "hxx")
region_extract_minimap(uuids, "env", "hxx")
region_extract_minimap(uuids, "gag", "hxx")
region_extract_minimap(uuids, "vif", "hxx")
region_extract_minimap(uuids, "nef", "hxx")
region_extract_minimap(uuids, "vpu", "hxx")
region_extract_minimap(uuids, "vpr", "hxx")
region_extract_minimap(uuids, "tat_exon1", "hxx")
region_extract_minimap(uuids, "tat_exon2", "hxx")
region_extract_minimap(uuids, "rev_exon1", "hxx")
region_extract_minimap(uuids, "rev_exon2", "hxx")

##merge exons from rev and tat (special case due to two codons)
exon_merge(uuids, "rev")
exon_merge(uuids, "tat")

# ---------------------------
# 4. Codon alignment with MACSE
# ---------------------------

for(region in c("env", "pol", "gag", "vif", "nef", "vpu", "vpr")){

  uuids = read.csv(paste0("blastmeta/blastmeta_",region,".csv")) %>%
    filter(!is.na(Alignment.Length)) %>%
    select(QueryID)
  if(nrow(uuids) == 0){
    print(paste0("No sequences for ", region))
    next
  }
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

  if(nrow(uuids) == 0){
    print(paste0("No sequences for ", region))
    next
  }

  foreach(x = unique(uuids$uuid)) %dopar% {

    codon_align(x, region)

  }
}

print(sessionInfo())