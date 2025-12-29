#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=3.6G        #physical memory limit
#$ -l h_vmem=3.6G      #virtual memory limit
#$ -l h_rt=100:00:00   #Wall time
#$ -j y                #merges STDOUT and STDERR
#$ -N base_recal1      #job name
#$ -t 1-27             #Number of files
#$ -tc 15              #maximal number of files running
#$ -pe smp 4
#$ -R y


export JAVA_HOME=/share/apps/jdk1.8.0_131
export PATH=$JAVA_HOME/bin:$PATH

#Load Picard tools
export PATH=/share/apps/genomics/picard-2.20.3/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/picard-2.20.3/bin:${LD_LIBRARY_PATH}
#Load SAM tools
export PATH=/share/apps/genomics/samtools-1.9/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/bin:${LD_LIBRARY_PATH}
#Load JAVA
export PATH=/share/apps/jdk1.8.0_131/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/jdk1.8.0_131/bin:${LD_LIBRARY_PATH}
#Load GATK
export PATH=/share/apps/genomics/gatk-4.1.2.0:${PATH}
export PATH=/share/apps/genomics/gatk-4.1.2.0:${LD_LIBRARY_PATH}
#Load pytho
export PATH=/share/apps/python-3.6.9/lib:${PATH}
export LD_LIBRARY_PATH=/share/apps/python-3.6.9/lib:${LD_LIBRARY_PATH}

#add any additional modules or paths here
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:${PATH}
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}


#https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR
cd /your/input/directory/

inputfile=( $(ls *_sorted.bam) ) 
k1="${inputfile[$SGE_TASK_ID-1]}"

genome=/directory/to/genome/GRCh38.d1.vd1.fa
SNPs=/directory/to/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
INDELS=/directory/to/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


java "-Xmx2G" -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar AddOrReplaceReadGroups \
        I=./$k1 O=./recal/${k1%%_sorted.bam}"_RG.bam" \
        ID=./${k1%%_sorted.bam} PU=B LB=DNA PL=ILLUMINA SM=./${k1%%_sorted.bam} 


#continue with script 2
        
