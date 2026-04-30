FROM ubuntu:24.04@sha256:84e77dee7d1bc93fb029a45e3c6cb9d8aa4831ccfcc7103d36e876938d28895b

RUN apt-get update && apt-get upgrade -y

# Install everything in root home directory, run everything as root.
WORKDIR /root

RUN apt-get install -y wget 

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py39_25.9.1-3-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /root/miniconda
ENV PATH="/root/miniconda/bin:$PATH"

# 'RUN' all subsequent shell commands in this Dockerfile with bash
SHELL [ "/bin/bash", "-c" ]
ARG CONDA_ALWAYS_YES=true

# Add conda channels
RUN conda config --prepend channels conda-forge
RUN conda config --prepend channels bioconda
RUN conda config --prepend channels r

ENV MAMBA_LOG_LEVEL=debug

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# R packages
RUN conda install -y -n base --solver=libmamba \
    -c conda-forge -c bioconda -c r \
    r-base=4.4 \
    bioconductor-biostrings=2.74.0 \
    r-doParallel=1.0.17 \
    r-iterators=1.0.14 \
    r-foreach=1.5.2 \
    r-tidyselect=1.2.1 \
    r-magrittr=2.0.3 \
    r-cli=3.6.5 \
    r-glue=1.8.0 \
    r-pillar=1.11.0 \
    r-vctrs=0.7.1 \
    r-codetools=0.2-20 \
    r-lifecycle=1.0.5 \
    r-data.table=1.17.8 \
    r-rlang=1.1.7 \
    r-fastDummies=1.7.5 \
    r-Rcpp=1.1.0 \
    r-Matrix \
    r-qdapTools=1.3.7 \
    r-svMisc=1.4.3 \
    r-glmnet=4.1-10 \
    r-dplyr=1.2.0 \
    r-rBLAST=1.3.1 \
    r-devtools=2.4.6 \
    r-tidyr=1.3.2 \
    r-seqinr=4.2-36 \
    bioconductor-decipher=3.2.0 \
    r-tidyverse=2.0.0 \
    r-aer=1.2.16 \
    r-ape=5.8 \
    r-phangorn=2.12.1 \
    r-ggrepel=0.9.6 \
    r-phytools=2.5.2 \
    r-phylolm=2.6.5 \
    r-lme4=1.1.38 \
    r-tidytable=0.11.1 \
    r-lmerTest=3.2-0 \
    r-brglm2=1.0.1 \
    r-lamw==2.2.7 \
    && conda clean -afy

# Some of the R packages we need are not available/have conflicts in conda, so we install them from source.
RUN R -e "\
  tryCatch( \
    install.packages('https://cran.r-project.org/src/contrib/ramcmc_0.1.2.tar.gz', repos = NULL, type = 'source'), \
    error=function(e) \
      install.packages('https://cran.r-project.org/src/contrib/Archive/ramcmc/ramcmc_0.1.2.tar.gz', repos = NULL, type = 'source') \
  )"

RUN R -e "\
  tryCatch( \
    install.packages('https://cran.r-project.org/src/contrib/adaptMCMC_1.5.tar.gz', repos = NULL, type='source'), \
    error=function(e) \
      install.packages('https://cran.r-project.org/src/contrib/Archive/adaptMCMC/adaptMCMC_1.5.tar.gz', repos = NULL, type = 'source') \
  )"

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/POUMM/POUMM_2.1.8.tar.gz", repos = NULL, type = "source")'

# Non R tools
RUN conda install -y -n base --solver=libmamba \
    -c conda-forge -c bioconda -c r \
    macse==2.07 \
    mafft==7.525 \
    plink==1.90b7.7 \
    plink2==2.0.0a.6.9 \
    eigensoft==8.0.0 \
    iqtree==3.1.1 \
    minimap2==2.30 \
    && conda clean -afy