###################
# Input folders:
# Note: sample name must be ffffffff-ffff-ffff-ffff-ffffffffffff
# - raw_bam // contains sample.bam
# - hxb2    // contains reference.fa(!sta)
# - map_fasta_cons  // contains consensus sequence of sample.bam, reference consensus for bam file, .fa

# Output folders (must already exist):
# - sorted_index_bam
# - ref_reads_to_hypermut
# - reads_to_hypermut
# - bamsummary
# - fastq_filtered
# - final_filtered_bam
# - unmapped_filtered_bam
###################

##Library
message("Loading Library................................................")
library(stringr)
library(tidyverse)
library(reticulate)
library(Rsamtools) # load explicitly for scanBamFlag helper function
##library(deepSNV) # bam2R - if loaded explicitly, seems to create problems: Error in (function (classes, fdef, mtable)  :   unable to find an inherited method for function 'select' for ignature '"data.frame"' Calls: %>% -> select -> <Anonymous>

##Local Function Library
##insert deletions into matched string
inject <- function(string, index, replacement){
  index = index - order(index) + 1
  stringi::stri_sub_replace_all(string, from = index,
                                to = index-1,
                                replacement = replacement)
}

g2a_muts = function(seq,ref){
  seq = strsplit(seq,"")[[1]]
  ref = strsplit(ref,"")[[1]]
  sum(ref[seq == "A"] == "G")
}

gsref_muts = function(ref){
  ref = strsplit(ref,"")[[1]]
  sum(ref == "G")
}


## Main Script
message("Main Script................................................")


##datapath to samtools
path_to_sam = "samtools"
##datapath to hypermut python script
path_to_hyper = "./hypermut.py"
#########################

##sammtools doesnt work in R, shell issue...
##copy into terminal
##run to sort bam file
i="bam_file_name"

  system(paste0(path_to_sam," sort ./raw_bam/",i,".bam -o ./sorted_index_bam/",i,"_sorted.bam"))
  system(paste0(path_to_sam," index ./sorted_index_bam/",i,"_sorted.bam"))
  

  ##filter for unmapped reads
  param <- Rsamtools::ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")

  ##exclude unmapped reads and save it in destination
  bam_bam <- Rsamtools::filterBam(paste0("./sorted_index_bam/",i,"_sorted.bam"), destination = paste0("./unmapped_filtered_bam/",i,"_sorted.bam"), param=param)

  ##read bam file with excluded unmapped reads
  bam_terst = Rsamtools::scanBam(paste0("./unmapped_filtered_bam/",i,"_sorted.bam"))
  
  initial_formating = as.data.frame(bam_terst[[1]][["seq"]])
  initial_formating$startpos = bam_terst[[1]][["pos"]]
  initial_formating$width = bam_terst[[1]][["qwidth"]]
  initial_formating$cigar = bam_terst[[1]][["cigar"]]
  
  initial_formating = initial_formating %>% mutate(chars = str_count(cigar,"[A-Z]+"))
  
  initial_formating = initial_formating %>% mutate(counts = paste0(str_extract_all(cigar, "[0-9]+"))) %>% 
    mutate(counts = gsub("\"|)","",counts)) %>%
    mutate(counts = gsub("c\\(","",counts)) %>%
    separate(counts, into = paste("count",1:max(.$chars), sep = "_"), sep = ",") %>%
    mutate(charcs = paste0(str_extract_all(cigar, "[A-Z]+"))) %>% 
    mutate(charcs = gsub("\"|)","",charcs)) %>%
    mutate(charcs = gsub("c\\(","",charcs)) %>%
    separate(charcs, into = paste("charsc",1:max(.$chars), sep = "_"), sep = ",") 
  
  initial_formating = initial_formating %>%
    mutate_at(vars(grep("count",names(.))),~ifelse(is.na(.),0,.)) %>%
    mutate_at(vars(grep("charsc",names(.))),~ifelse(is.na(.),"",.)) %>%
    mutate(id = rownames(.))
  
  fullcigar = initial_formating %>%
    select(id,grep("charsc",names(.)),grep("count",names(.)))
  
  cigarcount = (length(fullcigar)-1)/2
  for(y in 2:(cigarcount+1)) {
    fullcigar = cbind(fullcigar,str_dup(fullcigar[,y],fullcigar[,y+cigarcount]))
    colnames(fullcigar)[length(fullcigar)] = paste0("cigar",y-1)
  }
  
  fullcigar = fullcigar %>%
    unite(fullcigar,grep("cigar",names(.)), sep = "") %>%
    mutate(fullcigar = gsub(" ","",fullcigar)) %>%
    select(id,fullcigar)
  
  initial_formating = initial_formating %>% merge(.,fullcigar,by = "id", all = TRUE)
  
  ##takes forever
  ##inject deleitons as gaps into the read sequence
  deletions_insert = initial_formating %>%
    rowwise() %>%
    mutate(haha = inject(x,as.data.frame(gregexpr(pattern ='D',fullcigar))[,1],"-"))
  
  ##takes forever
  ##only extract matches/missmatches and deletions from the read sequence
  alignment_ready = deletions_insert %>%
    rowwise() %>%
    mutate(haha = paste0(str_sub(haha, as.data.frame(gregexpr(pattern ='M|D',fullcigar))[,1],as.data.frame(gregexpr(pattern ='M|D',fullcigar))[,1]),collapse = ""))
  
  
##quit(status = 0) # early exit just to see

  ##reference consensus for bam file
  fasta_msa = as.data.frame(readDNAStringSet(paste0("./map_fasta_cons/",i,".fa")))
  
  ##original
  fasta_org = fasta_msa %>%
    mutate(n = nchar(x)) %>%
    separate(x, into = paste("pos", 0:(max(.$n)), sep = ""), sep = "") %>%
    select(-pos0,-n) %>% mutate(id_org = rownames(.)) %>% pivot_longer(!id_org,names_to = "position_org", values_to = "base_org")
  
  ##read
  foralign = c(readDNAStringSet(paste0("./map_fasta_cons/",i,".fa")),readDNAStringSet(paste0("./hxb2.fa")))
  aligned = muscle::muscle(foralign)

  ##aligned with hxb2
  fasta_align = as.data.frame(as(aligned, "DNAStringSet")) %>%
    mutate(n = nchar(x)) %>%
    separate(x, into = paste("pos", 0:(max(.$n)), sep = ""), sep = "") %>%
    select(-pos0,-n) %>% mutate(id = rownames(.)) %>% pivot_longer(!id,names_to = "position", values_to = "base")
  
  ##take seq of interest and combine with original position
  fasta_test = fasta_align %>% filter(id != "K03455.1") %>% filter(base != "-")
  fasta_align_org_align = cbind(fasta_org,fasta_test)
  ##take hxb2 and merge with just above
  fasta_hxb2 = fasta_align %>% filter(id == "K03455.1")
  colnames(fasta_hxb2) = c("idhxb2","hxb2pos","hxb2base")
  
  fasta_align_three = merge(fasta_align_org_align,fasta_hxb2, by.x = "position", by.y = "hxb2pos", all = TRUE) %>%
    filter(!is.na(base)) %>% mutate(position = as.numeric(gsub("pos","",position))) %>% arrange(position)
  
  reference_fa = paste0(fasta_align_three$hxb2base,collapse = "")
  
  alignment_ready = alignment_ready %>%
    rowwise() %>%
    mutate(hahahaha = paste0(str_sub(reference_fa,startpos,startpos+nchar(haha)-1),collapse = ""))
  
  alignment_ready = alignment_ready %>%
    rowwise() %>%
    mutate(g2as = g2a_muts(haha,hahahaha)) %>%
    mutate(gs = gsref_muts(hahahaha)) %>%
    mutate(mutgs = g2as/gs)
  
  hypermut_ready = alignment_ready %>%
    rowwise() %>%
    mutate(midpos = startpos + round((nchar(haha)/2),0))
  
  
  ##calculate hypermutations as in hypermut from los alamos
  ##write read sequence and reference sequence in file
  write.table(paste(">",rownames(hypermut_ready),"\n",hypermut_ready$haha,"\n\n",collapse = '', sep = ""),
              paste0("./reads_to_hypermut/",i,"_ex_reads.fa"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(paste(">",rownames(hypermut_ready),"\n",hypermut_ready$hahahaha,"\n\n",collapse = '', sep = ""),
              paste0("./ref_reads_to_hypermut/",i,"_ex_ref_reads.fa"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  ##run python for hypermut algorithm 
  # note: we save the result of this python script to a csvso that we see the intermediate results
  hyyperr = as.data.frame(system(paste0(
    #"conda run -n hypermut ",
    "python ", # this also works
    path_to_hyper," ./reads_to_hypermut/",
                                        i,"_ex_reads.fa ./ref_reads_to_hypermut/",
                                        i,"_ex_ref_reads.fa | tee ./reads_to_hypermut/hyyperr.stdout.csv"), intern = TRUE))  # expected columns: (num_matched_muts,num_potential_muts,num_matched_ctrls,num_potential_ctrls,oddsratio,p)

  # for troubleshooting hyyperr.csv
  write.csv(hyyperr, "./reads_to_hypermut/hyyperr.as.data.frame.csv")

  colnames(hyyperr) = "muts"
  
  hyyperr = hyyperr %>% separate(col = muts,sep = ",", into = paste0("V",seq(1:6))) 
  ##read hypermut results
  ##extract hypermut p value
  hypermut_ready$hyper_p = as.numeric(gsub(")","",hyyperr$V6))
  hypermut_ready$bam = i
  
  ##save bam summary
  write.csv(hypermut_ready,paste0("./bamsummary/",i,"rna_ref_bam_summary.csv"))
  dyn_thr = read.csv("./dynamic_threshold.1.csv", header = TRUE) # expected columns: thr,pos
  hypermut_filter = hypermut_ready %>% merge(.,dyn_thr[,c("thr","pos")], by.x = "midpos", by.y = "pos", all.x = TRUE, all.y = FALSE) %>% 
    mutate(threshold = ifelse(is.na(thr),0.8,thr)) %>%
    mutate(filterflag = hyper_p > threshold) %>% group_by(id) %>% arrange(as.numeric(id))

  ##exclude unmapped reads and save it in destination
  BF = BamFile(
    paste0("./unmapped_filtered_bam/",i,"_sorted.bam"),
    yieldSize = 2000000 # TODO this should be dynamically parametrizeable via options.json ... there might be more number of max reads eventually...
  )
  path_to_final_bam = paste0("./final_filtered_bam/",i,"_final.bam")
  bam_bam <- Rsamtools::filterBam(BF,
                                  destination = path_to_final_bam, 
                                  param=ScanBamParam(what = scanBamWhat()),filter = hypermut_filter$filterflag)
  
  # Shows how man reads where left:
  countBam(bam_bam)

  # convert .bam to .fastq
  system(paste0(path_to_sam, " fastq ", path_to_final_bam, " > ./fastq_filtered/",i,"_filtered.fastq"))

  # save bam_frequencies.csv
  bam_frequencies <- deepSNV::bam2R(
    file = path_to_final_bam, 
    chr = "unknown_unknown_unknown_K03455.1:0:0:0", start = 1, stop=10000, q = 25
  )
  write.csv(bam_frequencies, "./final_filtered_bam/bam_frequencies.csv") 

  print(paste0("success ",i))
