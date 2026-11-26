#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=3.6G  
#$ -l h_vmem=3.6G   
#$ -l h_rt=200:00:00  
#$ -j y   
#$ -N make_bam1   
#$ -t 1-27    
#$ -tc 11  
#$ -pe smp 8
#$ -R y

#Load BWA
export PATH=/share/apps/genomics/bwa-0.7.17:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/bwa-0.7.17:${LD_LIBRARY_PATH}
#Load Picard tools
export PATH=/share/apps/genomics/picard-2.20.3/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/picard-2.20.3/bin:${LD_LIBRARY_PATH}
#Load SAM tools
export PATH=/share/apps/genomics/samtools-1.9/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/bin:${LD_LIBRARY_PATH}
#Load JAVA
export PATH=/share/apps/jdk1.8.0_131/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/jdk1.8.0_131/bin:${LD_LIBRARY_PATH}
#Load python
export PATH=/share/apps/python-3.6.9/lib:${PATH}
export LD_LIBRARY_PATH=/share/apps/python-3.6.9/lib:${LD_LIBRARY_PATH}

#add any additional modules or paths here
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:${PATH}
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}


#the genome file here should be indexed using bwa
genome=/directory/to/genome/GRCh38.d1.vd1.fa

#cd into the directory where the trimmed files are
cd /your/input/directory/trim


inputfile=( $( ls *_R1_001_val_1.fq.gz ) )
#This assigns the first item in the array as element zero (note the -1). Couting starts at 0
x1="${inputfile[$SGE_TASK_ID-1]}"
#Assign the mathcing R2 element to x2
x2=${x1%%_R1_001_val_1.fq.gz}"_R2_001_val_2.fq.gz"

#unzip both of these veriables
gunzip $x1
gunzip $x2

#assign the unzipped files to new variables.
y1=${x1%%.gz}
y2=${x2%%.gz}

#sample name, used for all output files
sample=${sample}

#Align reads to the genome using bwa
#mem is used for reads ? 70bps in length
#Remember to make an output directory
bwa mem -M -t 7 $genome $y1 $y2 > /your/directory/${sample}".sam" 
#Convert sam to bam
samtools view --threads 8 -bS /your/directory/${sample}".sam" -o /your/directory/${sample}".bam" 


#Create an uBAM file. This is an unmapped bam file
java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar FastqToSam \
        F1=$y1 F2=$y2 \
        TMP_DIR=/your/directory/temp \
        O=/your/directory/${sample}"_unmapped.bam" \
        SM=$sample

rm /your/directory/${sample}".sam" 

java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar MergeBamAlignment \
        TMP_DIR=/your/directory/temp \
        ALIGNED=/your/directory/${sample}".bam"  \
        UNMAPPED=/your/directory/${sample}"_unmapped.bam" \
        O=/your/directory/${sample}"_merged.bam" \
        R=$genome

###Remove unneccsary files          
rm /your/directory/${sample}".bam"
rm /your/directory/${sample}"_unmapped.bam"       

java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar MarkDuplicates \
        TMP_DIR=/your/directory/temp \
        I=/your/directory/${sample}"_merged.bam" \
        O=/your/directory/${sample}"_duplicates.bam" \
        M=/your/directory/${sample}"_duplicates_metrics.txt"

    
java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar SortSam \
        TMP_DIR=/your/directory/temp \
        I=/your/directory/${sample}"_duplicates.bam"\
        O=/your/directory/${sample}"_sorted.bam" \
        SORT_ORDER=coordinate  
        
rm /your/directory/${sample}"_merged.bam" 
rm /your/directory/${sample}"_duplicates.bam"



