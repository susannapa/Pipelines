#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=3.6G
#$ -l h_vmem=3.6G
#$ -l h_rt=100:00:00
#$ -j y
#$ -N gatk_working
#$ -t 1-27 #Number of files
#$ -tc 10
#$ -pe smp 4
#$ -R y

# Load JAVA
export JAVA_HOME=/share/apps/jdk1.8.0_131
export PATH=/share/apps/jdk1.8.0_131/bin:${PATH}
# Load Picard tools
export PATH=/share/apps/genomics/picard-2.20.3/bin:${PATH}
java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar -h
# Load SAM tools
export PATH=/share/apps/genomics/samtools-1.9/bin:${PATH}
# Load GATK
export PATH=/share/apps/genomics/gatk-4.1.2.0:${PATH}
# Load python
source /share/apps/source_files/python/python-3.6.9.source


cd /your/input/directory/

inputfile=( $(ls *_sorted.bam) )
k1="${inputfile[$SGE_TASK_ID-1]}"

genome=/directory/to/genome/GRCh38.d1.vd1.fa
SNPs=/directory/to/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
INDELS=/directory/to/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

gatk BaseRecalibrator \
    -I ./recal/${k1%%_sorted.bam}"_RG.bam" \
        -R $genome \
        --known-sites $SNPs \
        --known-sites $INDELS \
        -O ./recal/${k1%%_sorted.bam}"_recal_data.table"

gatk ApplyBQSR \
        -R $genome \
        -I ./recal/${k1%%_sorted.bam}"_RG.bam" \
        --bqsr-recal-file ./recal/${k1%%_sorted.bam}"_recal_data.table" \
        -O ./recal/${k1%%_sorted.bam}"_recal.bam"
