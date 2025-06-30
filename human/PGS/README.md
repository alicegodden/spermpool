# Polygenic Risk Score (PRS) Analysis Pipeline

This repository outlines a pipeline for conducting Polygenic Risk Score (PRS) analysis, from initial BAM file processing and variant calling to PRS calculation and downstream analysis.

---

## Pipeline Steps

### Step 1: Data Preparation and LiftOver

This initial step involves processing BAM files to generate VCFs, performing a `LiftOver` from GRCh38 to GRCh37, and standardizing chromosome naming.

#### 1.1 Variant Calling from BAM Files

BAM files were generated with `mpileup` and `bcftools call`. The following command shows how to generate VCF files from your aligned BAMs:

```bash
bcftools mpileup -a AD,ADF,ADR -r "${ncbi_chr}" -q 30 -Q 30 -f /filepath/GRCh38_latest_genomic.fna "${bam_file}" | \
bcftools call -mv -f GQ -Oz -o "$output_file"

1.2 Performing LiftOver to GRCh37
Polygenic Risk Scores (PGS) are typically calculated against the GRCh37 reference genome. Since our BAM files were aligned to GRCh38, a LiftOver step is necessary to convert the VCFs to GRCh37 coordinates.

First, ensure your GRCh37 reference genome (e.g., hg19.fa) is indexed:
```bash
samtools faidx hg19.fa

And generate a sequence dictionary:
```bash
gatk CreateSequenceDictionary \
  -R hg19.fa \
  -O hg19.fa.dict

