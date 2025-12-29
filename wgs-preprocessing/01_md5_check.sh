#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=1.6G       # physical memory limit
#$ -l h_vmem=1.6G     # virtual memory limit
#$ -l h_rt=966:00:00  # wall time
#$ -j y               # merge STDOUT and STDERR
#$ -N md5sum          # job name
#$ -R y
#$ -t 1-1             # tot number of files
#$ -tc 1              # max number of jobs running at once

#add any modules or paths here
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:${PATH}
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

#folder containing FASTQ files
INPUT_DIR="./FASTQ"

cd $INPUT_DIR

md5sum *.fastq.gz > calculated_md5_fastq.txt

#then compare the calculated_md5_fastq.txt file with the provided MD5 checksums file.
