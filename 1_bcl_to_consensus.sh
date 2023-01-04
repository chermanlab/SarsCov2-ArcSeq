#!/usr/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=60GB
#SBATCH -p x86
#SBATCH --output=<insert output directory>
#SBATCH -e <insert error directory>

#Activate conda environment
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh



conda activate arqseq



#Convert bcl files to fastq. Edit file locations for reads and sample sheet.
BCL2FASTQ_CMD="/gpfs0/apps/x86/anaconda3/bin/bcl2fastq"
cd #<setwd to run folder>
srun $BCL2FASTQ_CMD -R ./ -o ./fq1 \
--use-bases-mask Y*,I6N*,Y*,Y* \
--minimum-trimmed-read-length 5 \
--mask-short-adapter-reads 0 \
--sample-sheet ./<sample sheet directory> \
--stats-dir ./Stats/ \
--reports-dir ./Reports/ \



cd #<setwd to fastq folder>
echo "Concatenating files..."
gunzip *.gz
#Concatenate the fastq files from different lanes. For the NextSeq550, there are 4 lanes or data. 
#Data should be organized into 3 reads based on bcl2fastq commands above.
cat alpha_covid19_S1_L001_R1_001.fastq alpha_covid19_S1_L002_R1_001.fastq alpha_covid19_S1_L003_R1_001.fastq alpha_covid19_S1_L004_R1_001.fastq > alpha_covid19_R1.fastq

cat alpha_covid19_S1_L001_R2_001.fastq alpha_covid19_S1_L002_R2_001.fastq alpha_covid19_S1_L003_R2_001.fastq alpha_covid19_S1_L004_R2_001.fastq > alpha_covid19_R2.fastq

cat alpha_covid19_S1_L001_R3_001.fastq alpha_covid19_S1_L002_R3_001.fastq alpha_covid19_S1_L003_R3_001.fastq alpha_covid19_S1_L004_R3_001.fastq > alpha_covid19_R3.fastq



echo "Extracting UMIs..."
#Extract the UMIs and append to read headers using umi_tools
umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNN --stdin=alpha_covid19_R3.fastq --read2-in=alpha_covid19_R1.fastq --stdout=UMI_alpha_covid19_R1.fastq --read2-stdout
umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNN -I alpha_covid19_R3.fastq -S UMI_alpha_covid19_R3.fastq
umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNN --stdin=alpha_covid19_R3.fastq --read2-in=alpha_covid19_R2.fastq --stdout=UMI_alpha_covid19_R2.fastq --read2-stdout
umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=UMI_alpha_covid19_R2.fastq --read2-in=UMI_alpha_covid19_R1.fastq --stdout=UMI2_alpha_covid19_R1.fastq --read2-stdout
umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=UMI_alpha_covid19_R2.fastq --read2-in=UMI_alpha_covid19_R3.fastq --stdout=UMI2_alpha_covid19_R3.fastq --read2-stdout



cd #<setwd to fastqs>
echo "Creating cDNA consensus uBAMs"
picard FastqToSam \
F1=UMI2_alpha_covid19_R1.fastq \
F2=UMI2_alpha_covid19_R3.fastq \
O=alpha_covid19_ubam.bam \
SM=alpha_covid19_ubam_21July2021



cp alpha_covid19_ubam.bam #<end directory>
cd #<end directory>
echo "Reformatting unaligned BAM files..."
# collate the bam files
srun samtools collate alpha_covid19_ubam.bam -o alpha_covid19_ubam_collate.sam

# convert first '_' to '#' w/ awk
srun samtools view -h alpha_covid19_ubam_collate.sam | awk '{ sub(/\_/, "#"); print }' > alpha_covid19_ubam_collate_awk.sam

# convert second '_' to '-' w/ awk
srun samtools view -h alpha_covid19_ubam_collate_awk.sam | awk '{ sub(/\_/, "-"); print }' > alpha_covid19_ubam_collate_awk2.sam

# append RX header w/ awk
srun samtools view -h alpha_covid19_ubam_collate_awk2.sam | awk '/^@/ {print;next} {N=split($1,n,"#");print $0 "\tRX:Z:" n[N]}' > alpha_covid19_ubam_collate_awk2_RX.sam

# revert to compressed .bam format
srun samtools view -h alpha_covid19_ubam_collate_awk2_RX.sam -o alpha_covid19_ubam_collate_awk2_RX.bam



echo "Making consensus calls..."
# run consensus caller
srun python UnifiedConsensusMaker_Duplex_ARC_seqmerge.py --input alpha_covid19_ubam_collate_awk2_RX.bam --prefix alpha_covid19 --numAlignmentReads 0 --write-PCRcs
#srun python UnifiedConsensusMaker_Duplex_ARC_seqmerge.py --input alpha_covid19_ubam_collate_awk2_RX.bam --prefix alpha_covid19_maxmem1 --numAlignmentReads 0 --maxmem1 1000000 --write-PCRcs
#srun python UnifiedConsensusMaker_Duplex_ARC_seqmerge.py --input alpha_covid19_ubam_collate_awk2_RX.bam --prefix alpha_covid19_cutoff1 --numAlignmentReads 0 --maxmem1 1000000 --cutoff1 .9 --write-PCRcs



echo "Finished"
srun sleep 60 

