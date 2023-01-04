#!/usr/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=50GB
#SBATCH -p weka
#SBATCH --output=./bwa2muts_Ecoli_%j.out
#SBATCH -e ./bwa2muts_Ecoli_%j.err



#Activate conda environment
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh



conda activate ARCseq_2.3



echo "Index the referene genome..."
bwa index -a is Ecoli_NC_000913.3.fasta
echo "Create .dict file for reference genome..."
samtools dict Ecoli_NC_000913.3.fasta > Ecoli_NC_000913.3.dict



#gunzip *cDNAcs.fq.gz
echo "Map the consensus reads..."
bwa mem -R '@RG\tID:220704_NB551658_0188_AHFFVCBGXM\tSM:Omicron_calu3_aln_Ecoli\tLB:tARCseq\tPL:ILLUMINA' Ecoli_NC_000913.3.fasta Omicron_sarscov2_calu3_read1_cDNAcs.fq Omicron_sarscov2_calu3_read2_cDNAcs.fq | samtools sort -o Omicron_sarscov2_calu3_aln_Ecoli.bam
samtools index Omicron_sarscov2_calu3_aln_Ecoli.bam
samtools idxstats Omicron_sarscov2_calu3_aln_Ecoli.bam

echo "Generate mapping stats..."
samtools flagstat Omicron_sarscov2_calu3_aln_Ecoli.bam > stats_Omicron_sarscov2_calu3_aln_Ecoli

echo "Filtering mapped reads..."
#filter out unmapped reads with -F4 option. I'm not sure why a reference fasta is needed here.
samtools view -h -F4 Omicron_sarscov2_calu3_aln_Ecoli.bam | samtools view -Sb -T Ecoli_NC_000913.3.fasta - > Omicron_sarscov2_calu3_aln_Ecoli.filt.bam
samtools index Omicron_sarscov2_calu3_aln_Ecoli.filt.bam

echo "Clipping and GATK alignment correction... "
#picard AddOrReplaceReadGroups INPUT=Ecoli.sort.filt.bam OUTPUT=Ecoli_readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

bcftools mpileup -f Ecoli_NC_000913.3.fasta Omicron_sarscov2_calu3_aln_Ecoli.filt.bam > Omicron_sarscov2_calu3_aln_Ecoli.vcf
#bcftools --help



conda deactivate

#gatk --help
#Indel realignment
srun java -Xmx8G -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R Ecoli_NC_000913.3.fasta -I Omicron_sarscov2_calu3_aln_Ecoli.filt.bam -known Omicron_sarscov2_calu3_aln_Ecoli.vcf -o Omicron_sarscov2_calu3_aln_Ecoli.intervals
srun java -Xmx8G -jar GenomeAnalysisTK.jar -T IndelRealigner -R Ecoli_NC_000913.3.fasta -I Omicron_sarscov2_calu3_aln_Ecoli.filt.bam -targetIntervals Omicron_sarscov2_calu3_aln_Ecoli.intervals -known Omicron_sarscov2_calu3_aln_Ecoli.vcf -o Omicron_sarscov2_calu3_aln_Ecoli.realign.bam
srun java -Xmx8G -jar GenomeAnalysisTK.jar -T ClipReads -I Omicron_sarscov2_calu3_aln_Ecoli.realign.bam -R Ecoli_NC_000913.3.fasta --cyclesToTrim "1-6,146-150" --clipRepresentation SOFTCLIP_BASES -o Omicron_sarscov2_calu3_aln_Ecoli.realign_clipped.bam


conda activate arqseq


#Filter our low quality reads <q37. Output .bam with header.
samtools view -bh -q37 Omicron_sarscov2_calu3_aln_Ecoli.realign_clipped.bam > Omicron_sarscov2_calu3_aln_Ecoli.realign_clipped_q37.bam
samtools index Omicron_sarscov2_calu3_aln_Ecoli.realign_clipped_q37.bam


echo "Creating mpileup variant file..."
#Call variants.
#samtools mpileup -O -s --reference yeast_S288C.fasta sort_WTyeast_cDNAcs.bam  -o muts_WTyeast
samtools mpileup -B -d 1000000 -f Ecoli_NC_000913.3.fasta Omicron_sarscov2_calu3_aln_Ecoli.realign_clipped_q37.bam > muts_Omicron_sarscov2_calu3_aln_Ecoli.pileup
#bcftools mpileup -B -d 1000000 -f Ecoli_NC_000913.3.fasta WT_sarscov2_aln_Ecoli.realign_clipped_q37.bam > muts_WT_sarscov2_aln_Ecoli_bcftools.vcf


conda deactivate

echo "Tabulating variants..."
cat muts_Omicron_sarscov2_calu3_aln_Ecoli.pileup | python2 count-muts_KB.py -d 3 -C 0.5 > muts_Omicron_sarscov2_calu3_aln_Ecoli



#python --help
conda activate ARCseq_R_env

#R --version

Rscript Ecoli_covid_Rscript.R

echo "Finished"

srun sleep 60 
