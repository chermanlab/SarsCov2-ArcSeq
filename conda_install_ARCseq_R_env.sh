#!/usr/bin/sh
#SBATCH -o install.out
#SBATCH -e install.err
#SBATCH --job-name="Conda install"
#SBATCH -p x86

source /home/BCM/ccbradle/conda_bgs_14oct20/mini3/etc/profile.d/conda.sh

conda create -n ARCseq_R_env
conda activate ARCseq_R_env
conda install -c bioconda -c conda-forge \
    r-base=4.1.0 \
    r-Hmisc=4.6 \
    r-tidyverse=1.3.1 \
    r-xfun=0.27

