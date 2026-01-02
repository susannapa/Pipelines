#! /bin/bash
#$ -S /bin/bash
#$ -l tmem=3.6G  
#$ -l h_vmem=3.6G  
#$ -l h_rt=24:00:00  
#$ -j y    
#$ -N fastqc_trim   
#$ -t 1-27  # tot number of files (As we have forward and reverse reads it is the total number of samples=number.files/2)
#$ -tc 10  # max number of jobs running at once (don't flood the cluster!)

#export modules and path for the following:
#FastQC
export PATH=/share/apps/genomics/FastQC-0.11.8:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/FastQC-0.11.8:${LD_LIBRARY_PATH}
#Java
export PATH=/share/apps/openjdk-12.0.2/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/openjdk-12.0.2/bin:${LD_LIBRARY_PATH}
#cutadapt
export PATH=/share/apps/genomics/cutadapt-2.5/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/cutadapt-2.5/bin:${LD_LIBRARY_PATH}
#folder for python library
export LD_LIBRARY_PATH=/share/apps/python-3.6.4-shared/lib:${LD_LIBRARY_PATH}
#trim galore
export PATH=/share/apps/genomics/TrimGalore-0.6.3:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/TrimGalore-0.6.3:${LD_LIBRARY_PATH}

cd /your/input/directory/FASTQ

#lists the first reads of each sample as an array
inputfile=( $( ls *_R1_001.fastq.gz  ) )

#This assigns the first item in the array as element zero (note the -1). Couting starts at 0
x1="${inputfile[$SGE_TASK_ID-1]}"

#Assign the mathcing R2 element to x2
x2=${x1%%_R1_001.fastq.gz}"_R2_001.fastq.gz"

#first fastqc: 
    #t: number of threads to run
    #o: output directory
        fastqc -t 1 -o /your/output/directory/qc $x1  
        fastqc -t 1 -o /your/output/directory/qc $x2

 #then trimgalore: 
    #-j: number of cores
    #-q: trim low quality ends from read with a qulity score less than 20
    #--length: discards read that are less than 100nts in length after trimming
    #--trim-n: removes ns from eaitehr side of read
    #--path_to_cutadapt: specify a path to the Cutadapt executabl
    #--fastqc: run fastqc in default mode on teh trimmed files
    #-o where the trimed data will be outputted to. I suggest you create a trim folder within the Run1/FASTQ folder
    #-paired: the R1 and R2 files

        trim_galore -j 1 -q 20 --length 100 --trim-n \
        --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt \
        --fastqc \
        -o /your/output/directory/trim \
        --paired $x1 $x2


#then inspect the generated HTML reports. For each sample, we need to inspect the fastqc.html file from the qc and the trim folder, so after QCs and after trimming
