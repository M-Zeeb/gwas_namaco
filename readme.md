# GWAS NAMACO

## Repository for dockerized pipelines of HIV GWAS and heriatibility analyses for the manuscript published in xxx...

### Performs HIV GWAS and Heritability analysis 
#### Requirements (see folder input for example (created with "create_example_data.R")):
 1. A folder with frequencey files from smaltalign output

 2. A phenotype file

 3. File with outcomes defined

 4. Files with covariables defined

### If docker is installed, complete example analysis can be run by executing "run.sh"
### Without docker manual installation of all software is required (see DOCKERFILE)

### Issues/TODO

 1. Dynamic parallelisation in "foreach-steps" depending on available cores

 2. Mamba package installations often get stuck at random points...
