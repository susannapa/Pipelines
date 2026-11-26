#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
WGS preprocessing: FASTQ -> trimmed reads -> BAM -> BQSR -> gVCF
same steps as the bash scripts 02-06, one sample per task
Author: Susanna Pagni
*/

params.reads  = "/your/input/directory/FASTQ/*_R{1,2}_001.fastq.gz"
params.genome = "/path/to/GRCh38.d1.vd1.fa"      // indexed with bwa, samtools faidx and a .dict
params.dbsnp  = "/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
params.mills  = "/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.outdir = "results"


process fastqc_trimgalore {
    tag "$sample"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.{html,zip,txt}"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_val_{1,2}.fq.gz"), emit: trimmed
    path "*.{html,zip,txt}"

    script:
    """
    fastqc -t 1 ${reads}
    trim_galore -j 1 -q 20 --length 100 --trim-n --fastqc --paired ${reads}
    """

    stub:
    """
    touch ${sample}_R1_001_val_1.fq.gz ${sample}_R2_001_val_2.fq.gz ${sample}_fastqc.html
    """
}

process make_bam {
    tag "$sample"
    cpus 8

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_sorted.bam")

    script:
    """
    #align, adding the read group here (instead of AddOrReplaceReadGroups later)
    bwa mem -M -t ${task.cpus} -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:DNA" \
        ${params.genome} ${reads} | samtools sort -@ ${task.cpus} -o ${sample}_aligned.bam

    picard MarkDuplicates I=${sample}_aligned.bam O=${sample}_sorted.bam M=${sample}_duplicates_metrics.txt
    rm ${sample}_aligned.bam
    """

    stub:
    """
    touch ${sample}_sorted.bam
    """
}

process base_recalibration {
    tag "$sample"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_recal.bam")

    script:
    """
    gatk BaseRecalibrator -I ${bam} -R ${params.genome} \
        --known-sites ${params.dbsnp} --known-sites ${params.mills} \
        -O ${sample}_recal_data.table

    gatk ApplyBQSR -I ${bam} -R ${params.genome} \
        --bqsr-recal-file ${sample}_recal_data.table -O ${sample}_recal.bam
    """

    stub:
    """
    touch ${sample}_recal.bam
    """
}

process call_gvcf {
    tag "$sample"
    publishDir "${params.outdir}/gvcf", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    path "${sample}.g.vcf.gz"

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.genome} -I ${bam} \
        -O ${sample}.g.vcf.gz -ERC GVCF
    """

    stub:
    """
    touch ${sample}.g.vcf.gz
    """
}


workflow {
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

    fastqc_trimgalore(reads)
    make_bam(fastqc_trimgalore.out.trimmed)
    base_recalibration(make_bam.out)
    call_gvcf(base_recalibration.out)
}
