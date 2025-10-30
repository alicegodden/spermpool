# 🎯 TE & Variant Analysis  
*(human / GRCh38)*

## Overview  
This sub-directory contains scripts for transposable element (TE) insertion analysis, variant calling, phenogram visualizations, and related data-filtering workflows tailored to the human genome (GRCh38).  

---

## Directory Contents  

| File | Description |
|------|-------------|
| [`phenogram.py`](./phenogram.py) | Generates chromosome-level phenograms incorporating TE insertions, gene labels, centromere locations, and chromatin accessibility tracks. |
| [`filter_retroseq.py`](./filter_retroseq.py) | Filters output from a TE-calling pipeline (e.g., RetroSeq) by read depth, confidence metrics, and experimental vs control group status. |
| [`class_chart.py`](./class_chart.py) | Produces summary charts of TE insertions by class/family (e.g., LINE, SINE, LTR) based on filtered data. |
| [`TE_Chrom_count.py`](./TE_Chrom_count.py) | Counts TE insertion events per chromosome and visualizes their distribution to detect hotspots and genome-wide patterns. |
| [`vcfinput_phenogram.py`](./vcfinput_phenogram.py) | Creates phenogram visualizations from VCF files (e.g., annotated variant calls from SnpEff/SnpSift workflows) — useful for overlaying variant impact on chromosome layout. |
| [`filter.py`](./filter.py) | Filters annotated VCF files by variant impact level or specific functional terms (e.g., “High” impact, GO-term hits) for downstream visualization. |

---

## Getting Started  

### Prerequisites  
- **Python** ≥ 3.11  
- Required Python libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`, `argparse`, `biopython`  
- Input files must match expected formats:  
  - VCF files with annotation (for `vcfinput_phenogram.py`, `filter.py`)  
  - RetroSeq output (for `filter_retroseq.py`, `phenogram.py`)  
  - Reference supporting files: e.g., centromere positions (`chrcen.txt`), gene coordinates (`genes.txt`), chromosome-end positions (`chrom_end.txt`)  

### Example Workflow  
```bash
# 1. Filter TE calls from RetroSeq output
python filter_retroseq.py
Filtering Retroseq outputs for quality

# 2. Count insertions per chromosome
python TE_Chrom_count.py
Requires: chrom_end.txt and centromere_grch38

# 3. Generate phenogram visualization
python phenogram.py
Needs access to these files:
    --genes genes.txt \
    --chrcen chrcen.txt \
    --chrom_ends chrom_end.txt \
    --output phenogram_plot.png

# 4. For variant-impact visualization:
python filter.py 
python vcfinput_phenogram.py 
