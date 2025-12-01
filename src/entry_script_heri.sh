#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
trap 'echo FATAL ERROR EXIT CODE $? AT $0:$LINENO' ERR

workingdirectory=""
fasta_files=""
phenotype_files=""
ref_seq=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --wd)
      workingdirectory="$2"
      shift 2
      ;;
    --f)
      fasta_files="$2"
      shift 2
      ;;
    --r)
      ref_seq="$2"
      shift 2
      ;;
    --p)
      phenotype_files="$2"
      shift 2
      ;;
    --*) # catch unknown options
      echo "Unknown option: $1"
      exit 1
      ;;
    *) # catch positional args
      echo "Unexpected argument: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$fasta_files" ]]; then
  #echo "Error: --s for raw ngs file source is required"
  #exit 1
  echo "No fasta file source given, default used"
  fasta_files="heri_input/fasta_files"
fi

if [[ -z "$ref_seq" ]]; then
  #echo "Error: --s for reference sequence is required"
  #exit 1
  echo "No reference seq given, default used"
  ref_seq="heri_input/HIV1_REF_2023_pol_DNA.fa"
fi

if [[ -z "$phenotype_files" ]]; then
  #echo "Error: --p for phenotype/covariate file is required"
  #exit 1
  echo "No phenotype/covariate file given, default used"
  phenotype_files="heri_input/phenotype_files"
fi

if [[ -z "$workingdirectory" ]]; then
  echo "No working directory given, default directory used"
  workingdirectory="heri_output"
fi

echo "workingdirectory = $workingdirectory"
echo "fasta_files = $fasta_files"
echo "phenotype_files = $phenotype_files"

echo "############################"
echo "Successfully defined basic folders"
echo "Creating subfolders and copy input files"
echo "############################"


#Create folder structure
mkdir -p ${workingdirectory}/{phenotype_files,results,fasta_files,for_alignment_nt,alignment_nt,references_seqs,phylos,poumm,phylo_clusters}
mkdir -p ${workingdirectory}/poumm/{bs_summaries_pmm,bs_summaries_poumm}

# copy raw files, phenotype data, and reference sequence files in working directory
cp $fasta_files/*  ${workingdirectory}/fasta_files/
cp $phenotype_files/* ${workingdirectory}/phenotype_files/
cp $ref_seq  ${workingdirectory}/references_seqs/ref_seq.fa

# absolute path to source files
src_path=$(realpath "src/")

#Change directory to working directory
cd ${workingdirectory}

# Phenotypes
mapfile phenos < ./phenotype_files/phenotype_outcome.txt

echo "############################"
echo "Successfully created all subfolders and input files copied"
echo "Creating alignments"
echo "############################"

cat ./fasta_files/*.fa > ./for_alignment_nt/for_alignment_whole_temp1.fa
awk '/^>/{flag=(index($0,"HXB2")>0)?1:0} !flag' ./for_alignment_nt/for_alignment_whole_temp1.fa > ./for_alignment_nt/for_alignment_temp2.fa
sed '/^>/! s/[-*!]//g' ./for_alignment_nt/for_alignment_temp2.fa > ./for_alignment_nt/for_align.fa
rm ./for_alignment_nt/*temp*.fa

echo pwd

mafft --localpair --maxiterate 1000  --thread 64 \
    --add ./for_alignment_nt/for_align.fa \
    ./references_seqs/ref_seq.fa \
    >  ./alignment_nt/aligned.fa

cp ./alignment_nt/aligned.fa ./phylos/

echo "############################"
echo "Successfully created alignment"
echo "Creating tree"
echo "############################"

iqtree -s ./phylos/aligned.fa \
    -alrt 1000 \
    -B 1000 \
    -T AUTO \
    -wbtl

echo "############################"
echo "Successfully created tree"
echo "Calculating heritability MIXED"
echo "############################"

for pheno in "${phenos[@]}"; do
  Rscript ${src_path}/heritability_mixed.R $pheno
done

echo "############################"
echo "Successfully calculated heritability mixed"
echo "Preparing heritability data for poumm"
echo "############################"

for pheno in "${phenos[@]}"; do
  Rscript ${src_path}/prepare_poumm.R $pheno
done

echo "############################"
echo "Successfully created Heritability data"
echo "Calculating Heritability POUMM"
echo "############################"

for pheno in "${phenos[@]}"; do
  Rscript ${src_path}/poumm.R $pheno
done

echo "############################"
echo "Successfully calculated Heritability POUMM"
echo "Creating figures"
echo "############################"

Rscript ${src_path}/heritability_visualization.R

echo "############################"
echo "Everything Done!!"
echo "############################"
