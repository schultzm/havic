conda create -n havic_env python=3.8 R
export PATH=/bin:/usr/bin
echo "channels:
  - bioconda
  - conda-forge
  - r
  - defaults
channel_priority: strict
auto_activate_base: false" > ~/.condarc
echo "install.packages('BiocManager', repos='https://cran.ms.unimelb.edu.au/')
install.packages('pheatmap',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('qualpalr',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('ape',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('caper',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('diversitree',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('geiger',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('nlme',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('OUwie',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('phangorn',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('phytools',repos='https://cloud.r-project.org',quiet=TRUE)
install.packages('igraph',repos='https://cloud.r-project.org',quiet=TRUE)
BiocManager::install('Rsamtools')
BiocManager::install('RCurl')
BiocManager::install('GenomicAlignments')
BiocManager::install('Biostrings')
BiocManager::install('Rsamtools')
BiocManager::install('ggtree')
BiocManager::install('Biobase')
BiocManager::install('matrixStats')
BiocManager::install('tidyverse')" > rinstall.R
R --file=rinstall.R
conda install samtools minimap2 iqtree clusterpicker graphviz
conda update --all
havic depcheck
# iqtree                      : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/iqtree)
# minimap2                    : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/minimap2)
# R                           : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/R)
# samtools                    : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/samtools)
havic test


