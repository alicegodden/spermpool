
# Polygenic Risk Score (PGS) Analysis Pipeline

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
```

#### 1.2 Performing LiftOver to GRCh37
Polygenic Risk Scores (PGS) are typically calculated against the GRCh37 reference genome. Since our BAM files were aligned to GRCh38, a LiftOver step is necessary to convert the VCFs to GRCh37 coordinates.

First, ensure your GRCh37 reference genome (e.g., hg19.fa) is indexed:

```bash
samtools faidx hg19.fa
```

And generate a sequence dictionary:

```bash
gatk CreateSequenceDictionary \
  -R hg19.fa \
  -O hg19.fa.dict
```

Then, perform the LiftOver using GATK's LiftoverVcf:

```bash
gatk LiftoverVcf \
  -I Sample_GRCh38_ensembl_renamed.vcf.gz \
  -O GRCh37_D1T0.vcf \
  -R Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa \ # Reference: [https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz) (v2015-11-27)
  -CHAIN GRCh38_to_GRCh37.chain.gz \ # Obtained from: [https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/)
  -REJECT rejected.vcf
```

#### 1.3 Renaming Chromosomes
The LiftOver output VCF files may use RefSeq chromosome identifiers (e.g., "NC_"). The following sed script renames these to Ensembl-style (e.g., "1", "X", "Y"):

```bash

bcftools view input.vcf.gz | \
sed -e 's/^NC_000001.11/1/' \
    -e 's/^NC_000002.12/2/' \
    -e 's/^NC_000003.12/3/' \
    -e 's/^NC_000004.12/4/' \
    -e 's/^NC_000005.10/5/' \
    -e 's/^NC_000006.12/6/' \
    -e 's/^NC_000007.14/7/' \
    -e 's/^NC_000008.11/8/' \
    -e 's/^NC_000009.12/9/' \
    -e 's/^NC_000010.11/10/' \
    -e 's/^NC_000011.10/11/' \
    -e 's/^NC_000012.12/12/' \
    -e 's/^NC_000013.11/13/' \
    -e 's/^NC_000014.9/14/' \
    -e 's/^NC_000015.10/15/' \
    -e 's/^NC_000016.10/16/' \
    -e 's/^NC_000017.11/17/' \
    -e 's/^NC_000018.10/18/' \
    -e 's/^NC_000019.10/19/' \
    -e 's/^NC_000020.11/20/' \
    -e 's/^NC_000021.9/21/' \
    -e 's/^NC_000022.11/22/' \
    -e 's/^NC_000023.11/X/' \
    -e 's/^NC_000024.10/Y/' | \
bgzip > output.chromnames.vcf.gz
```

### Step 2: Extracting Allele Frequencies
This step extracts Alternate Allele Frequencies (AAF) from your VCF files. This is crucial for subsequent PRS calculations.

The following bash script uses bcftools and awk to extract chromosome, position, reference allele, alternate allele, and the calculated AAF into a TSV file.

```bash

#!/bin/bash

vcf_file="GRCh37_sample.vcf"
tsv_file="GRCh37_sample.aaf.tsv" # Alternate Allele Frequency output

if [ ! -f "$vcf_file" ]; then ## Check if the VCF file exists
    echo "ERROR: VCF file not found: $vcf_file" >&2
    exit 1
fi

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' "$vcf_file" | \ ## Use bcftools to extract the relevant information, handling potential errors
awk '
    BEGIN {
        FS = "\t";
        OFS = "\t";
        print "CHROM", "POS", "REF", "ALT", "AAF";
    }
    {
        if ($0 ~ /CHROM\tPOS\tREF\tALT\t\[AD\]/) next; # Skip the header line from bcftools

        # Check if AD field exists
        if (NF >= 5 && $5 != ".") {
            if ($5 ~ /^[0-9]+,[0-9]+$/) {
                split($5, ad, ",");
                total = ad[1] + ad[2];
                aaf = (total > 0) ? ad[2] / total : "NA";
                print $1, $2, $3, $4, aaf;
            } else {
                print "WARNING: AD field has unexpected format: " $5 " in line " NR > "/dev/stderr";
                print $1, $2, $3, $4, "NA"; # Output NA for AAF
            }
        } else {
            print "WARNING: AD field is missing in line " NR > "/dev/stderr";
            print $1, $2, $3, $4, "NA"; # Output NA for AAF
        }
    }
' > "$tsv_file"

echo "Finished processing VCF file and writing to $tsv_file"
```

### Step 3: Downloading and Processing PGS Files
This step involves downloading a substantial dataset of Polygenic Score (PGS) files from the Pan-UK Biobank and filtering them to retain relevant columns for analysis.

#### 3.1 Downloading PGS Files
PGS files were downloaded from the Pan-UK Biobank phenotype manifest, accessible here. The URLs for these files are stored in a file named urls.txt.

The following Slurm-compatible bash script handles the parallel download and initial filtering of these files:

```bash

#!/bin/bash
#SBATCH --job-name=download_ukb
#SBATCH --output=download_%a.out
#SBATCH --error=download_%a.err
#SBATCH -p compute
#SBATCH --array=1-7228%100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00

echo "Starting download and filtering for task ${SLURM_ARRAY_TASK_ID}..."

process_file() { # Function to download and filter a single file
    url="$1"
    filename=$(basename "$url")
    basename="${filename%.tsv.bgz}"
    logfile="${basename}.log" # Create a log file

    echo "Downloading $filename"
    # Remove -q for verbose output, add error handling
    wget "$url" -O "$filename" 2>>"$logfile" || {
        echo "ERROR: Download failed for $filename" >>"$logfile"
        return 1 # Return non-zero to indicate failure
    }

    echo "Filtering $filename"
    zcat "$filename" | awk '
        BEGIN {
            FS=OFS="\t";
            # Print header to a separate file for debugging (optional)
            # print > (FILENAME ".header")
        }
        NR==1 {
            # Map columns by exact name (case-insensitive)
            for (i=1; i<=NF; i++) {
                col = tolower($i)
                if (col == "chr") chr=i
                else if (col == "pos") pos=i
                else if (col == "ref") ref=i
                else if (col == "alt") alt=i
                else if (col == "af_meta_hq") af=i
                else if (col == "beta_meta_hq") beta=i
                else if (col == "se_meta_hq") se=i
                else if (col == "neglog10_pval_meta_hq") pval=i
            }

            # Fail early if any column is missing
            if (!(chr && pos && ref && alt && af && beta && se && pval)) {
                print "ERROR: Missing one or more required columns in " FILENAME > "/dev/stderr"
                return 1
            }
            # Print the desired header
            print "chr", "pos", "ref", "alt", "af_meta_hq", "beta_meta_hq", "se_meta_hq", "neglog10_pval_meta_hq"
            next
        }
        {
            print $chr, $pos, $ref, $alt, $af, $beta, $se, $pval
        }
    ' | bgzip > "${basename}.filtered.tsv.bgz" 2>>"$logfile" || {
        echo "ERROR: bgzip failed for ${basename}.filtered.tsv.bgz" >>"$logfile"
        return 1
    }

    # Check if the output file is empty
    if [ ! -s "${basename}.filtered.tsv.bgz" ]; then
        echo "ERROR: Output file ${basename}.filtered.tsv.bgz is empty" >>"$logfile"
        return 1
    fi

    echo "Cleaning up $filename"
    rm "$filename"
    return 0 # Explicitly return 0 for success
}

declare -a urls # Read the URLs from the file into an array.
while IFS= read -r url; do
    urls+=("$url")
done < "urls.txt"

if ((SLURM_ARRAY_TASK_ID > ${#urls[@]})); then # Check if the array index is valid
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is out of bounds (max is ${#urls[@]})" >&2
    exit 1
fi

if ! process_file "${urls[$((SLURM_ARRAY_TASK_ID - 1))]}"; then # Process the URL corresponding to the array task ID.
    echo "ERROR: process_file failed for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Finished processing task ${SLURM_ARRAY_TASK_ID}."
```

#### 3.2 Filtering Downloaded PGS Files
After downloading, the PGS .tsv.bgz files are further filtered to retain specific columns relevant for PRS calculation, including: chr, pos, ref, alt, beta_meta, beta_meta_hq, beta_EUR, and neglog10_pval_EUR.

```bash

zcat categorical-6141-both_sexes-8.tsv.bgz | awk -F'\t' -v OFS='\t' '
NR==1 {
  n=split("chr pos ref alt beta_meta beta_meta_hq beta_EUR neglog10_pval_EUR", want, " ");
  for (i=1; i<=NF; i++) col[$i]=i;
  for (i=1; i<=n; i++) if (want[i] in col) keep[++k]=col[want[i]];
  print_header="";
  for (i=1; i<=k; i++) print_header = print_header (i>1 ? OFS : "") want[i];
  print print_header;
}
NR>1 {
  line = "";
  for (i=1; i<=k; i++) line = line (i>1 ? OFS : "") $keep[i];
  print line;
}
' | gzip -c > output_filt.tsv.bgz

```


### Step 4: Matching VCF and PGS Data
This step involves matching the loci in your sample VCF files with the overlapping loci in the downloaded and filtered PGS .tsv.bgz files. This is a crucial step for preparing the data for PRS calculation.

The following Slurm-compatible bash script iterates through the filtered PGS files and executes a Python script to perform the matching.

```bash
#!/bin/bash
#SBATCH --job-name=match_aaf
#SBATCH --output=logs/match_%A_%a.out
#SBATCH --error=logs/match_%A_%a.err
#SBATCH -p compute
#SBATCH --array=0-7219%300
#SBATCH --mem=8G
#SBATCH --time=07:00:00

module load python/anaconda/2023.07/3.11.4

FILES=($(ls *.tsv.bgz)) # Get list of files
FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$FILE" .tsv.bgz)

python match_aaf_retain_all_tsv_columns.py "$FILE" "matched_parts/${BASENAME}_matched.tsv" # Run matching
```

The match_aaf_retain_all_tsv_columns.py script is responsible for this matching process. It multiplies the Alternate Allele Frequency (AAF) from your VCFs with the beta_EUR value from the PGS files to generate a risk score for each phenotype.

### Analysis in R & Python
This section outlines the Python scripts used for further analysis of the calculated polygenic risk scores.

phenocodes: Contains all the phenotype codes and their descriptions/associated files. This can be a text file or a small script to manage this information.

phenotype_score_difference_positive.R # Plotting traits with a PGS score above 0, showing the top 40 PGS scored traits

histogram.R, histogram_raw_score.R,  histogram_raw_score_MC.R # histogram plotting raw and difference of PGS with different input file handling

heatmap.R # Plotting heatmaps

dumbbell_Custom_chunks.R # Dumbbell plots

dumbbell_plot_composite_MC_%.R

dumbbell_plot_composite_SU_%.R

pearsons_heatmaps_SU.py / pearsons_heatmaps_MC.py # Pearson correlation coefficient analysis of traits
