
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Nextflow wrapper for WGS preprocessing
Author: Susanna Pagni
*/

params.fastq_dir = "/your/input/directory/FASTQ"
params.work_dir  = "/your/project/work"
params.genome    = "/path/to/GRCh38.fa"

workflow {
    md5_check()
    fastqc_trimgalore()
    make_bam()
    base_recalibration()
    call_variants()
}

/*
----------------------------------
Processes
----------------------------------
*/

process md5_check {

    tag "MD5"
    
    script:
    """
    bash 01_md5_check.sh
    """
}

process fastqc_trimgalore {

    tag "FastQC_TrimGalore"

    script:
    """
    bash 02_fastqc_trimgalore.sh
    """
}

process make_bam {

    tag "MakeBAM"

    script:
    """
    bash 03_make_bam.sh
    """
}

process base_recalibration {

    tag "BQSR"

    script:
    """
    bash 04_base_recal_1.sh
    bash 05_base_recal_2.sh
    """
}

process call_variants {

    tag "VariantCalling"

    script:
    """
    bash 06_call_gvcf.sh
    bash 07_call_vcf_single.sh
    """
}
