#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=2:00:00
#$ -l h_vmem=2G,tmem=2G

export PATH=/share/apps/genomics/bcftools-1.9/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/bcftools-1.9/bin:${LD_LIBRARY_PATH}
export PATH=/share/apps/genomics/htslib-1.10/bin:${PATH}
export LD_LIBRARY_PATH=/share/apps/genomics/htslib-1.10/bin:${LD_LIBRARY_PATH}


cd /your/directory

# Read input file into array
readarray -t FILE < /path/to/your/files.txt
#files.txt looks like this:
#00340208189.vcf.gz
#00340224566.vcf.gz
#00340224610.vcf.gz

y="${FILE[$SGE_TASK_ID-1]}"

for index in ${!FILE[*]}; do
		bcftools view ${FILE[$index]}_merge_hf.vcf -f PASS $y > /path/${FILE[$index]}_merge_hf_filtered.vcf
		bgzip -c ${FILE[$index]}_merge_hf_filtered.vcf > ${FILE[$index]}_merge_hf_filtered.vcf.gz
    tabix -p vcf ${FILE[$index]}_merge_hf_filtered.vcf.gz
done

