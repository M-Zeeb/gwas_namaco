FROM ubuntu:20.04@sha256:f8f658407c35733471596f25fdb4ed748b80e545ab57e84efbdb1dbbb01bd70e

RUN echo "Dummy echo to force complete rebuild: 2023-08-11"

# TODO update puts some of the system in an unreproducible state... but without it, the urls apt knows might be outdated and the content be just gone, especially for fast changing dependencies (security relevant ones...)
RUN apt-get update && apt-get upgrade -y

# Install everything in root home directory, run everything as root.
WORKDIR /root

# Dependency of conda installation step:
# note: throughout this, we are assuming apt-get installed dependencies are very stable, no need to pin versions...
RUN apt-get install -y wget 

# Install conda environment manager, miniconda distribution. 
#
# see
# https://repo.continuum.io/miniconda
# alternatively: https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh -O miniconda.sh

RUN bash miniconda.sh -b -p /root/miniconda
ENV PATH="/root/miniconda/bin:$PATH"
RUN echo $PATH
RUN which conda
RUN conda --version
RUN conda list


# 'RUN' all subsequent shell commands in this Dockerfile with bash, to support >() syntax etc.
SHELL [ "/bin/bash", "-c" ]

# avoid 3 sec delay on conda/mamba install with "Confirm changes: [Y/n] "
ARG CONDA_ALWAYS_YES=true

# 'conda' always failed to install/resolve conflicts while installing smalt,
# so we choose to use the alternative 'mamba' package/environment manager instead
# alternatively try 0.24.0 https://github.com/mamba-org/mamba/issues/1728
RUN conda install mamba -n base -c conda-forge
RUN mamba config set extract_threads 1
# HACKFIX https://github.com/mamba-org/mamba/issues/1826#issuecomment-1196636463 ImportError: libarchive.so.13: cannot open shared object file: No such file or directory
# strangely, this does nothing and just reports > # All requested packages already installed.
RUN conda install libarchive -n base -c conda-forge

RUN which mamba
# RUN ldd "$(which mamba)" # not a dynamic executable
RUN mamba --version
RUN mamba list

# First add channels
RUN conda config --prepend channels conda-forge
RUN conda config --prepend channels bioconda
RUN conda config --prepend channels r

ENV MAMBA_LOG_LEVEL=debug

# hardcoded from known good build https://gitlab.shcs.ch/INF_NGS_DATABASE/medvir-smaltalign-container/-/jobs/170393
RUN mamba install r-base==4.4 -vv
RUN mamba install bioconductor-biostrings -vv
RUN mamba install r-tidytable==0.11
RUN mamba install r-doParallel==1.0.17
RUN mamba install r-iterators==1.0.14 
RUN mamba install r-foreach==1.5.2
RUN mamba install r-tidyselect==1.2.1
RUN mamba install r-magrittr==2.0.3
RUN mamba install r-cli==3.6.5
RUN mamba install r-glue==1.8.0 
RUN mamba install r-pillar==1.11.0
RUN mamba install r-vctrs==0.6.5
RUN mamba install r-codetools==0.2-20
RUN mamba install r-lifecycle==1.0.4
RUN mamba install r-data.table==1.17.8
RUN mamba install r-rlang==1.1.6
RUN mamba install r-fastDummies==1.7.5
RUN mamba install r-Rcpp==1.1.0 
RUN mamba install r-Matrix==1.7-3
RUN mamba install r-qdapTools==1.3.7
RUN mamba install r-svMisc==1.4.3
RUN mamba install r-glmnet==4.1-10
RUN mamba install r-dplyr
RUN mamba install r-rBLAST
RUN mamba install r-devtools -vv
RUN mamba install r-tidyr -vv
RUN mamba install r-seqinr -vv
RUN mamba install bioconductor-DECIPHER -vv
RUN mamba install r-tidyverse -vv
RUN mamba install r-AER -vv
RUN mamba install r-ape -vv
RUN mamba install r-phangorn -vv
RUN mamba install r-ggrepel -vv
RUN mamba install r-phytools -vv
RUN mamba install r-phylolm -vv
RUN mamba install r-lme4 -vv
RUN R -e 'install.packages("lmerTest", repos = "https://cloud.r-project.org")'
RUN R -e 'install.packages("POUMM", repos = "https://cloud.r-project.org")'

RUN mamba install bioconda::macse -vv

RUN mamba install mafft -vv

RUN mamba install bioconda::plink -vv
RUN mamba install bioconda::plink2 -vv

RUN mamba install eigensoft -vv

RUN mamba install bioconda::iqtree -vv