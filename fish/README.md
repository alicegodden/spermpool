# spermpool/fish

## 🚀 Overview  
Here are all the analysis and plotting scripts for Zebrafish (GRCz11/DanRer11) analysis

## Getting Started  
### Prerequisites  
- Python — version 3.11 or above
- R v4.1.1 or above

---

## 📁 Directory Contents

### 🎯 Variant Calling & Phenogram Visualization

#### `fish_phenogram.py`
Plots sperm pool variant calling data as a chromosome-level variant visualization.  
Useful for visualizing SNP/indel distributions and variant clusters across chromosomes.

**Inputs:**
- Processed VCF file(s)  
- Optional metadata (pool IDs, variant annotations)

**Outputs:**
- PNG/PDF phenogram plots of variant distribution  
- Optional overlay of SNP density or annotation categories

---

### 🧬 GO Term Visualization

#### `autobubble_goplot.py`
Automates the creation of **bubble plots** from GO term enrichment data (e.g., from *ShinyGO* CSV outputs).  
It standardizes visual style and scales for reproducibility.

**Inputs:**
- CSV file exported from [ShinyGO](http://bioinformatics.sdstate.edu/go/)

**Outputs:**
- Bubble plot summarizing GO term enrichment (Size = gene count, Color = p-value)

---

### 🧫 RetroSeq & Transposable Element (TE) Analysis

#### `testing_glm.R`
Performs **Generalized Linear Model (GLM)** analysis of **transposable element (TE)** counts filtered for:
- Read count thresholds  
- Confidence levels  
- Control group filtering (removes TEs present in control samples)

**Purpose:**
Used to identify TEs significantly enriched or depleted in experimental sperm pools.

---

#### `phenogram.py`
Generates a **comprehensive phenogram** integrating:
- RetroSeq-derived variant calls (VCF)
- Chromatin accessibility
- Gene annotation
- Centromere locations

**Inputs (DanRer11-compatible):**
- `Your.vcf` – RetroSeq output  
- `chrcen.txt` – Chromosome centromere positions  
- `genes.txt` – Gene coordinates and labels  
- `chrom_end.txt` – Chromosome end positions  

**Outputs:**
- Chromosome phenogram with genes, TE locations, and accessibility tracks

---

#### Supporting Scripts for RetroSeq Phenogram

| File | Purpose |
|------|----------|
| `filter_retroseq.py` | Filters RetroSeq output for confidence, read depth, and experimental relevance |
| `class_chart.py` | Generates summary charts by TE class/family |
| `TE_Chrom_count.py` | Counts and visualizes TE insertions per chromosome |

---

### 🧩 SnpEff / SnpSift Variant Filtering

#### `vcfinput_phenogram.py`
Plots VCF-derived variants (from SnpEff or SnpSift pipelines) as a **phenogram**, highlighting functional categories and chromosomal clustering.

**Inputs:**
- Annotated VCF file from SnpEff/SnpSift

**Outputs:**
- Chromosome plot highlighting variant locations by impact category

---

#### `filter.py`
Filters annotated VCF files by **variant impact** or **functional term**, e.g.:
- High impact variants only
- Variants with certain scores




