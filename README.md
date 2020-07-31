# HAVIC

[![CircleCI](https://circleci.com/gh/schultzm/HAVIC.svg?style=svg&circle-token=9d17418bb752aa29e07f95b09af106aef7cc6b02)](https://app.circleci.com/pipelines/github/schultzm/HAVIC)

Detect **H**epatitis **A** **V**irus **I**nfection **C**lusters from HAVNET amplicon sequences.  

## Usage

    havic
    havic depcheck
    havic test
    havic detect
    havic version



## Installation

Installation requires [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

    # install miniconda - see https://docs.conda.io/en/latest/miniconda.html
    export PATH=/bin:/usr/bin
    . $HOME/miniconda3/etc/profile.d/conda.sh
    conda init
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict --set auto_activate_base false
    conda create --yes --name havic_env python=3.8 R minimap2 iqtree clusterpicker graphviz
    conda activate havic_env
    conda update --yes --all #maybe not necessary
    conda install --yes -c bioconda samtools # This is installed separately to avoid the problems described at https://github.com/conda/conda/issues/8103
    echo "install.packages('BiocManager', repos='https://cran.ms.unimelb.edu.au/')
    BiocManager::install('GenomicAlignments')
    BiocManager::install('Biostrings')
    BiocManager::install('Rsamtools')
    BiocManager::install('ggtree')
    BiocManager::install('tidyverse')
    install.packages('Rcpp', repos='https://cloud.r-project.org') # see discussion here if \"make: 77 error\" occurs: https://github.com/cole-trapnell-lab/monocle3/issues/318#issuecomment-639521739
    install.packages('phytools',repos='https://cloud.r-project.org', type='source') " | R --no-save
    pip install git+https://github.com/schultzm/HAVIC.git
    havic depcheck
    havic test
