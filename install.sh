#!/bin/bash

# install miniconda - see https://docs.conda.io/en/latest/miniconda.html
conda init bash
source ~/.bashrc
export PATH=${PATH}:/bin:/usr/bin:/sbin:/usr/sbin
conda update conda
conda env create -f environment.yml
conda activate havic_env
havic version
conda list -n havic_env
arr=(hav_amplicon hav_wgs hav_pmc measles_wgs hiv_amplicon)
for i in ${arr[@]}
  do
    havic test ${i}
  done
