#!/bin/bash

# install miniconda - see https://docs.conda.io/en/latest/miniconda.html
conda init bash
source ~/.bashrc
git pull
export PATH=${PATH}:/bin:/usr/bin:/sbin:/usr/sbin
conda clean --all
conda update conda
conda install -c conda-forge mamba
mamba create -n havic_env python==3.9.0 r-base==4.0.3
conda activate havic_env
mamba env update --file environment.yml --prune
havic version
arr=(hav_amplicon hav_wgs hav_pmc measles_wgs hiv_amplicon)
for i in ${arr[@]}
  do
    echo havic test ${i}
  done
