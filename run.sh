#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
trap 'echo FATAL ERROR EXIT CODE $? AT $0:$LINENO' ERR

sudo docker build --progress plain --file Dockerfile -t shcsgwas:latest .

# Run the GWAS pipeline
sudo docker run \
    --volume "./output:/root/output" \
    --volume "./input:/root/input" \
    --volume "./src:/root/src" \
    --volume "./utility:/root/utility" \
    shcsgwas:latest /bin/bash /root/src/entry_script.sh

# Create fasta file input for heritability analysis
mkdir -p heri_input/fasta_files
mkdir -p heri_input/phenotype_files
cp input/phenotype_files/* heri_input/phenotype_files/

# Copy pol fasta sequences and filter out HXB2
cat output/codon_align/codon_align_pol_NT/* > heri_input/fasta_files/fasta_file_tmp.fa
awk 'BEGIN{RS=">"; ORS=""} NR>1 && tolower($1) !~ /hxb2/ \
    {print ">"$0}' heri_input/fasta_files/fasta_file_tmp.fa > \
    heri_input/fasta_files/filtered.fasta
# Rename fasta header from base_uuids to sample IDs
awk '
FNR==NR{
  a[$1]=$2
  next
}
($2 in a) && /^>/{
  print ">"a[$2]
  next
}
1
' heri_input/phenotype_files/pheno.names.txt FS="[> ]"  heri_input/fasta_files/filtered.fasta >\
    heri_input/fasta_files/fasta_file.fa

rm heri_input/fasta_files/filtered.fasta heri_input/fasta_files/fasta_file_tmp.fa

# Run the heritability pipeline
sudo docker run \
    --volume "./heri_output:/root/heri_output" \
    --volume "./heri_input:/root/heri_input" \
    --volume "./src:/root/src" \
    --volume "./utility:/root/utility" \
    shcsgwas:latest /bin/bash /root/src/entry_script_heri.sh
    