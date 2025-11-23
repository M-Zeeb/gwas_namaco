#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
trap 'echo FATAL ERROR EXIT CODE $? AT $0:$LINENO' ERR

workingdirectory=""
raw_files=""
phenotype_files=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --wd)
      workingdirectory="$2"
      shift 2
      ;;
    --s)
      raw_files="$2"
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

if [[ -z "$raw_files" ]]; then
  #echo "Error: --s for raw ngs file source is required"
  #exit 1
  echo "No raw ngs file source given, default used"
  raw_files="input/raw_seqs"
fi

if [[ -z "$phenotype_files" ]]; then
  #echo "Error: --p for phenotype/covariate file is required"
  #exit 1
  echo "No phenotype/covariate file given, default used"
  phenotype_files="input/phenotype_files"
fi

if [[ -z "$workingdirectory" ]]; then
  echo "No working directory given, default directory used"
  workingdirectory="output"
fi

echo "workingdirectory = $workingdirectory"
echo "raw_files = $raw_files"
echo "phenotype_files = $phenotype_files"

echo "############################"
echo "Successfully defined basic folders"
echo "Creating subfolders and copy input files"
echo "############################"

#Path to phenotypes

#ID col

#Outcome

#Covariables


#Create folder structure
mkdir -p ${workingdirectory}/{raw_files,phenotype_files,freq_files,fasta_files,blasted_seqs,blastmeta,codon_align,results,for_alignment_aa,for_alignment_nt,alignment_aa,alignment_nt,alignment_basefreq,pca,gwas_data}

# copy raw files, phenotype data, and reference sequence files in working directory
cp $raw_files/*  ${workingdirectory}/raw_files/
cp $phenotype_files/*  ${workingdirectory}/phenotype_files/
cp -a utility/references_seqs ${workingdirectory}/

# absolute path to source files
src_path=$(realpath "src/")

#Change directory to working directory
cd ${workingdirectory}

genes=("gag" "env" "pol" "rev" "tat" "vif" "vpr" "vpu" "nef")
exons=("rev_exon1" "rev_exon2" "tat_exon1" "tat_exon2")
# Loop over arrays to create subdirectories
#genes
for gene in "${genes[@]}"; do
    mkdir -p "blasted_seqs/blasted_seqs_$gene"
    mkdir -p "codon_align/codon_align_${gene}_NT"
    mkdir -p "codon_align/codon_align_${gene}_AA"
    mkdir -p "pca/pca_${gene}"
done
#exons
for exon in "${exons[@]}"; do
    mkdir -p "blasted_seqs/blasted_seqs_$exon"
done
#whole
    mkdir -p "pca/pca_whole"

echo "############################"
echo "Successfully created all subfolders and input files copied"
echo "Creating frequency files"
echo "############################"

# Sequencing formatting, NGS assembly (smalt align ouput) to structure frequency file, to consensus fasta file
Rscript ${src_path}/formatting_seqs.R "raw_files/" "freq_files/" "fasta_files/"

echo "############################"
echo "Successfully created frequency files"
echo "Performing blast, codon alignments, and APD score calculation"
echo "############################"

# Extract regions from each sequence with Blast and generate codon alignments 
Rscript ${src_path}/codon_align.R "freq_files/" "fasta_files/"

echo "############################"
echo "Successfully performed blast and codon alignments"
echo "Calculating APD scores"
echo "############################"

# Calculate apd for pol, gag, and env 
Rscript ${src_path}/apd_score.R "freq_files/" "fasta_files/"

echo "############################"
echo "Successfully calculated APD scores"
echo "Creating gene wise alignments"
echo "############################"

# Alignment Genes
for gene in "${genes[@]}"; do
    cat ./codon_align/codon_align_${gene}_AA/*.fa > ./for_alignment_aa/for_alignment_${gene}_temp1.fa
    awk '/^>/{flag=(index($0,"HXB2")>0)?1:0} !flag' ./for_alignment_aa/for_alignment_${gene}_temp1.fa > ./for_alignment_aa/for_alignment_${gene}_temp2.fa
    sed '/^>/! s/[-*!]//g' ./for_alignment_aa/for_alignment_${gene}_temp2.fa > ./for_alignment_aa/for_alignment_${gene}.fa
    rm ./for_alignment_aa/*temp*.fa
done

for gene in "${genes[@]}"; do
    mafft --localpair --maxiterate 1000 --thread 64 \
    --add ./for_alignment_aa/for_alignment_${gene}.fa \
        ./references_seqs/hxb2_${gene}_aa.fa \
    > ./alignment_aa/alignment_${gene}.fa
done

for gene in "${genes[@]}"; do
    Rscript ${src_path}/aa_to_nt.R $gene
done

# Alignment whole genome

echo "############################"
echo "Successfully created gene wise alignments"
echo "Creating whole genome alignment"
echo "############################"

cat ./fasta_files/*.fa > ./for_alignment_nt/for_alignment_whole_temp1.fa
awk '/^>/{flag=(index($0,"HXB2")>0)?1:0} !flag' ./for_alignment_nt/for_alignment_whole_temp1.fa > ./for_alignment_nt/for_alignment_whole_temp2.fa
sed '/^>/! s/[-*!]//g' ./for_alignment_nt/for_alignment_whole_temp2.fa > ./for_alignment_nt/for_alignment_whole.fa
rm ./for_alignment_nt/*temp*.fa

#mafft --localpair --maxiterate 1000 --thread 64 \
mafft --thread 32 --threadtb 8 \
    --add ./for_alignment_nt/for_alignment_whole.fa \
    ./references_seqs/hxb2_wga.fa \
    > ./alignment_nt/alignment_whole.fa

# PCA  
echo "############################"
echo "Successfully created whole genome alignments"
echo "Performing PCA"
echo "############################"

# Plink preformatting
for gene in "${genes[@]}" whole; do

    Rscript ${src_path}/plink.R $gene

    cd pca/pca_$gene
  
    ##PLINK
    plink2 --vcf vcf_${gene}.vcf \
    --keep-allele-order --allow-extra-chr --chr-set -1 --make-pgen --sort-vars \
    --out gwasp2
  
    plink2 --pfile gwasp2 \
                --pheno dummy_pheno_${gene}.txt \
                --nonfounders --allow-extra-chr --make-pgen \
                --out hivp2_pheno1
  
    ##remove with missing phenotypes
    plink2 --pfile hivp2_pheno1 \
                --require-pheno dummy_pheno --allow-extra-chr --make-pgen \
                --out hivp2_pheno2
  
     ##linkage correlation matrix hopefully
    plink2 --pfile hivp2_pheno2 \
                --geno 0.2 --nonfounders --make-bed \
                --out ${gene}_forcor
  
    ##linkage
    plink2 --pfile hivp2_pheno2 \
                --geno 0.2 --nonfounders --indep-pairwise 50 10 0.5 \
                --out gcta_${gene}_exclude
  
    ##save for gcta
    plink2 --pfile hivp2_pheno2 \
                --geno 0.2 --nonfounders \
                --exclude gcta_${gene}_exclude.prune.out --make-bed \
                --out gcta_$gene
  
    ##further for eigensoft
    plink2 --pfile hivp2_pheno2 \
                --mind 0.95 --nonfounders --allow-extra-chr  --make-pgen \
                --out hivp2_mind3
  
    plink2 --pfile hivp2_mind3 \
                --geno  0.99  --nonfounders --allow-extra-chr --make-pgen \
                --out hivp2_mind3
  
    ##make bed for eigensoft
    plink2 --pfile hivp2_mind3 \
                --nonfounders --allow-extra-chr --export ped \
                --out $gene

    #pheno type in column six as 1 (errors otherwise with "0" or "-9" values)
    awk 'BEGIN{OFS="\t"}{$6=1; print}' "${gene}.ped" > "${gene}.tmp" && mv "${gene}.tmp" "${gene}.ped"

  cd ../..

done

#PCA with smartpca

# activate environment for eigensoft
#source activate eigensoft

for gene in "${genes[@]}" whole; do

    cd pca/pca_$gene

    #Create parfile for covertf
    printf '%s\n' \
        "genotypename:    ${gene}.ped" \
        "snpname:         ${gene}.map" \
        "indivname:       ${gene}.ped" \
        "outputformat:    EIGENSTRAT" \
        "genotypeoutname: ${gene}.eigenstratgeno" \
        "snpoutname:      ${gene}.snp" \
        "indivoutname:    ${gene}.ind" \
        "familynames:     NO" > "${gene}_par.ped.eigenstrat"

    #Run convertf to create eigensoft compatible files
    convertf -p "${gene}_par.ped.eigenstrat"

    #Uniform population column for eigensoft
    awk 'BEGIN{OFS="\t"}{$3=1; print}' "${gene}.ind" > "${gene}.tmp" && mv "${gene}.tmp" "${gene}.ind"

    #Create parfile for smartpca
    printf '%s\n' \
        "genotypename:  ${gene}.eigenstratgeno" \
        "snpname:       ${gene}.snp" \
        "indivname:     ${gene}.ind" \
        "evecoutname:   ${gene}.evec" \
        "evaloutname:   ${gene}.eval" \
        "altnormstyle:  NO" \
        "numoutevec:    10" \
        "familynames:   NO" \
        "grmoutname:    grmjunk" > "${gene}_parfile"

    #Run smartpca
    smartpca -p "${gene}_parfile" > "${gene}_logfile_test_check"

    cd ../..

done

# activate environment for eigensoft
#conda deactivate

#GWAS 
for gene in "${genes[@]}"; do
    #pre formatting
    Rscript ${src_path}/gwas_preformatting.R $gene
    # run GWAS
    Rscript ${src_path}/gwas.R $gene
done

Rscript ${src_path}/gwas_visualization.R