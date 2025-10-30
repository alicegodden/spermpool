# Title: Human TE counts normalised by chromosome length
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def extract_chrom_pos_from_vcf(vcf_file):
    """Extract chromosome positions from the VCF file."""
    chrom_pos_dict = {}
    try:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue  # Skip header lines
                fields = line.strip().split('\t')
                if len(fields) < 2:  # Ensure there are enough fields
                    continue
                chromosome = fields[0]  # Correctly get chromosome
                position = int(fields[1])  # Get position as integer
                if chromosome not in chrom_pos_dict:
                    chrom_pos_dict[chromosome] = []
                chrom_pos_dict[chromosome].append(position)
    except Exception as e:
        print(f"Error reading {vcf_file}: {e}")
    return chrom_pos_dict

def read_chrom_lengths(chrom_file):
    """Read chromosome lengths from file and return a dictionary with lengths in Mb."""
    chrom_lengths = {}
    try:
        with open(chrom_file, 'r') as f:
            for line in f:
                chrom, length = line.strip().split()
                chrom_lengths[chrom] = int(length) / 1_000_000  # Convert to Mb
    except Exception as e:
        print(f"Error reading chromosome lengths from {chrom_file}: {e}")
    return chrom_lengths

# VCF files and sample names
vcf_files = [
    'raw_merged_donor_dedup_win100_gq200_fl8_chrom.vcf',
    'rcentral_merged_donor_dedup_minusRaw_win100_gq200_fl8_chrom.vcf',
    'outer_merged_donor_dedup_minusRaw_win100_gq200_fl8_chrom.vcf'
]

sample_names = ['Raw', 'Central', 'Outer']

# Read chromosome lengths
chrom_lengths = read_chrom_lengths('GRCh38_chrom_length.txt')

# Store normalized TE counts for each VCF file
all_te_counts = []

# Process each VCF file
for vcf_file in vcf_files:
    chrom_pos_dict = extract_chrom_pos_from_vcf(vcf_file)
    if chrom_pos_dict:  # Only process if valid data was extracted
        te_counts = {chrom: len(positions) / chrom_lengths.get(chrom, 1) for chrom, positions in chrom_pos_dict.items()}  # Normalize by chromosome length
        all_te_counts.append(te_counts)
    else:
        print(f"No valid data extracted from {vcf_file}.")

# Prepare data for plotting
chromosomes = list(chrom_lengths.keys())  # Ensure chromosome order is consistent
data = {sample_name: [te_counts.get(chrom, 0) for chrom in chromosomes] for sample_name, te_counts in zip(sample_names, all_te_counts)}

# Debugging output to check data structure
print("Data for plotting:")
for sample_name, te_counts in data.items():
    print(f"{sample_name}: {te_counts}")

# Color palette for the lines
rocket_palette = sns.color_palette("autumn", len(sample_names))

# Plotting
fig, ax = plt.subplots(figsize=(15, 8))

for sample_name, color in zip(sample_names, rocket_palette):
    if sample_name in data:
        ax.plot(chromosomes, data[sample_name], label=sample_name, color=color, marker='o')
    else:
        print(f"Warning: No data available for sample {sample_name}.")

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('TE counts per Mb', fontsize=18, fontweight='bold')
plt.title('MethylCellulose TEs per Chromosome (Normalized by Length) from Retroseq', fontsize=20, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.legend(title='Samples', fontsize=12, title_fontsize=14, loc='upper right')

# Save and display the plot
output_file = 'TE_counts_per_chromosome_per_MB_hapselDK_MC_retroseq.png'
plt.savefig(output_file, dpi=600)
plt.show()
