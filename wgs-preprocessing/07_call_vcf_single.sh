#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=3.6G   
#$ -l h_vmem=3.6G  
#$ -l h_rt=120:00:00  
#$ -j y    
#$ -N bamtovcf_single 
#$ -pe smp 4
#$ -t 1-4 
#$ -tc 4 

#this script works for single samples

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

#cd where the recal bams are located
cd /your/input/directory/recal

inputfile=( *recal.bam  )
k1="${inputfile[$SGE_TASK_ID-1]}"

genome=/directory/to/genome/GRCh38.d1.vd1.fa

gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R $genome \
     -I $k1 \
     -O /directory/vcf/${k1%%_.recal.bam}".vcf.gz" \
               
