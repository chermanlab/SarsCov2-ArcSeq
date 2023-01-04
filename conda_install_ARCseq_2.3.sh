#!/usr/bin/sh
#SBATCH -o install.out
#SBATCH -e install.err
#SBATCH --job-name="Conda install"
#SBATCH -p x86

source /home/BCM/ccbradle/conda_bgs_14oct20/mini3/etc/profile.d/conda.sh

conda create -n ARCseq_2.3
conda activate ARCseq_2.3
conda install -c bioconda -c conda-forge \
    bwa=0.7.17 \
    picard=2.18.29 \
    pysam=0.15.3 \
    python=3.6 \
    samtools=1.7 \
    gatk=3.8=7 \
    r=3.6 \
    bcftools \
    bedtools \
    umi_tools \
    openjdk=8.0.192


conda install -c dranew bcl2fastq

