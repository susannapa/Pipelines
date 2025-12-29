#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=3.8G          #physical memory limit
#$ -l h_vmem=3.8G        #virtual memory limit
#$ -l h_rt=2:00:00       #Wall time
#$ -j y                  #merges STDOUT and STDERR
#$ -N hardfil            #job name
#$ -pe smp 4
#$ -t 1-785              #total number of jobs you want to run
#$ -tc 50                #tot number of jobs running at the same time


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
#Load python
export PATH=/share/apps/python-3.6.9/lib:${PATH}
export LD_LIBRARY_PATH=/share/apps/python-3.6.9/lib:${LD_LIBRARY_PATH}
#load bcf tools
export PATH=/share/apps/genomics/bcftools-1.9/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/bcftools-1.9/bin:${LD_LIBRARY_PATH}

#add any additional modules or paths here
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:${PATH}
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}


# Read input file into array
readarray -t inputfile < /path/to/your/files.txt

#files.txt looks like this:
#00340208189.vcf.gz
#00340224566.vcf.gz
#00340224610.vcf.gz

# Read path file into array
readarray -t path < /path/to/your/file_path.txt

#file_path.txt looks like this:
#/path/to/00340208189.vcf.gz
#/path/to/00340224566.vcf.gz
#/path/to/00340224610.vcf.gz

y="${inputfile[$SGE_TASK_ID-1]}"
p="${path[$SGE_TASK_ID-1]}"

cd /your/input/directory/
genome=/directory/to/genome/GRCh38.d1.vd1.fa

#select SNPs
gatk --java-options "-Xmx4g" SelectVariants \
    -R $genome \
    -V $p \
    -select-type SNP \
    -O ${y%%.vcf.gz}"_SNP.vcf.gz"

#filter SNPs
gatk --java-options "-Xmx4g" VariantFiltration \
		-R $genome \
		-V /path/${y%%.vcf.gz}"_SNP.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O /path/${y%%.vcf.gz}"_SNP_hf.vcf.gz"


#select INDELS
gatk --java-options "-Xmx4g" SelectVariants \
    -R $genome \
    -V $p \
    -select-type INDEL \
    -O ${y%%.vcf.gz}"_indels.vcf.gz"

#filter INDELS
gatk --java-options "-Xmx4g" VariantFiltration \
      -R $genome \
	    -V /path/${y%%.vcf.gz}"_indels.vcf.gz" \
	    -filter "QD < 2.0" --filter-name "QD2" \
    	-filter "QUAL < 30.0" --filter-name "QUAL30" \
    	-filter "FS > 200.0" --filter-name "FS200" \
    	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"\
	    -O /path/${y%%.vcf.gz}"_indels_hf.vcf.gz"

#merge SNPs and indels
java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar MergeVcfs \
        I=/path/${y}"_SNP_hf.vcf.gz" \
        I=/path/${y}"_indels_hf.vcf.gz" \
        O=/path/${y}"_merge_hf.vcf"

