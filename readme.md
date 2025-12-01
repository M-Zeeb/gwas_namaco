# GWAS NAMACO

## Repository for dockerized pipelines of HIV GWAS and heriatibility analyses for the manuscript published in xxx...

### GWAS Requirements:
 1. A folder with frequencey files from smaltalign output

 2. A phenotype file

 3. File with outcomes defined

 4. FIles with covariables defined

### Heritability Requirements:
 1. A folder with consensus fasta files

 2. A phenotype file

 3. File with outcomes defined

 4. FIles with covariables defined

### TODO
 1. Generate/add example data

 2. Dynamic parallelisation in "foreach-steps" depending on available cores

 3. Mamba package installations often get stuck at random points...

 4. Generalize bam_hypermut_filter.R script

 5. Check necessary references
