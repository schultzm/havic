#!/bin/bash

# install miniconda - see https://docs.conda.io/en/latest/miniconda.html
conda init bash
source ~/.bashrc
export PATH=${PATH}:/bin:/usr/bin:/sbin:/usr/sbin
conda clean --all
conda update conda
conda create -n havic_env -c conda-forge python==3.9.0 r-base==4.0.3
conda activate havic_env
# conda env update --file environment.yml --prune
pip install -e .
havic version
conda list -n havic_env
arr=(hav_amplicon hav_wgs hav_pmc measles_wgs hiv_amplicon)
# for i in ${arr[@]}
#   do
    # echo havic test ${i}
#   done
