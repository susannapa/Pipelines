# Genomics Analysis Pipelines

This repository contains bioinformatics pipelines for processing WGS and variant calling data, designed for large-scale genomic studies. The scripts have been used in academic projects and are structured for reproducibility and HPC environments.

## Repository Structure

### wgs-preprocessing/

This folder contains scripts for pre-processing raw sequencing data (FASTQ) into BAM files and variant calls.

01_md5_check.sh – Validates file integrity using MD5 checksums
02_fastqc_trimgalore.sh – Performs quality control and trimming of raw FASTQ files
03_make_bam.sh – Aligns reads to reference genome and generates BAM files
04_base_recal_1.sh & 05_base_recal_2.sh – Performs base quality recalibration using GATK
06_call_gvcf.sh – Generates gVCF files per sample using HaplotypeCaller
07_call_vcf_single.sh – Calls VCFs from gVCFs for single samples

### vcf-hard-filtering/

This folder contains scripts to perform hard filtering of VCF files after variant calling.

01_vcf_hard_filter.sh – Selects SNPs and INDELs from raw VCFs; applies hard filters on quality metrics (QD, QUAL, FS, MQ, ReadPosRankSum, etc.) using GATK. Finally, it merges filtered SNP and INDEL files into a final VCF per sample.
02_filter_pass_variants.sh – Filters merged VCFs to include only variants with PASS status. Compresses and indexes the filtered VCFs with bgzip and tabix for downstream analyses.


## Usage

Modify paths in scripts to point to your input data, reference genome, and output directories.
Prepare sample lists as .txt files if using array jobs.
Submit scripts to your HPC scheduler (qsub / sbatch) or run locally with adjusted resource settings.

### Dependencies
GATK ≥4.1
bcftools ≥1.9
htslib ≥1.10
Picard Tools
BWA
Python3
Java 1.8+
Linux / HPC environment recommended
