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
#samtools                  1.4.1                         0    bioconda
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



#install.packages('pheatmap',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('qualpalr',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('ape',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('caper',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('diversitree',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('geiger',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('nlme',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('OUwie',repos='https://cloud.r-project.org',quiet=TRUE)
# install.packages('phangorn',repos='https://cloud.r-project.org',quiet=TRUE)#?this is in phytools
# install.packages('igraph',repos='https://cloud.r-project.org',quiet=TRUE)#? this is in phytools
# BiocManager::install('RCurl')

# BiocManager::install('Biobase')"
# #BiocManager::install('matrixStats')
# > rinstall.R
# R --file=rinstall.R
# conda update --all
# havic depcheck
# # iqtree                      : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/iqtree)
# # minimap2                    : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/minimap2)
# # R                           : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/R)
# # samtools                    : ok (/usr/local/Caskroom/miniconda/base/envs/havic_env/bin/samtools)
# havic test


