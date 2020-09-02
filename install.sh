# install miniconda - see https://docs.conda.io/en/latest/miniconda.html
export PATH=/bin:/usr/bin
. $HOME/miniconda3/etc/profile.d/conda.sh
conda env create -f environment.yml
conda init
source ~/.bashrc
conda activate havic_env
echo "install.packages('BiocManager', repos='https://cloud.r-project.org')
BiocManager::install('GenomicAlignments')
BiocManager::install('Biostrings')
BiocManager::install('Rsamtools')
BiocManager::install('ggtree')
BiocManager::install('tidyverse')
install.packages('pheatmap', repos='https://cloud.r-project.org')
install.packages('Rcpp', repos='https://cloud.r-project.org')
install.packages('phytools',repos='https://cloud.r-project.org', type='source') " | R --no-save
havic test
