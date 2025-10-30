# Title: Retroseq with donor
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd
import re
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

# Define file paths
vcf_input_file = 'use_T4_merged_donor_dedup_minusT0_win100_chrom_gq200_fl8.vcf'  # VCF file path
centromere_file = 'centromere_grch38.txt'  # Centromere positions file
genes_file = 'GRCh38_genes.txt'  # Genes information file
end_file = 'GRCh38_chrom_length.txt'  # Chromosome lengths file

# Read the input VCF file
data = pd.read_csv(vcf_input_file, sep='\t', comment='#', header=None)
data.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'T4files']

# Extract sample names from the last column
data['Sample'] = data['T4files'].apply(
    lambda x: re.search(r'source=([^\s]+)', x).group(1) if re.search(r'source=([^\s]+)', x) else None)

# Drop rows where Sample is NaN
data.dropna(subset=['Sample'], inplace=True)

# Function to map chromosomes
def get_chrom_number(chrom):
    if chrom == "chrX":
        return 23
    elif chrom == "chrY":
        return 24
    else:
        try:
            return int(chrom.replace("chr", ""))
        except ValueError:
            return None  # Exclude any other chromosomes

# Create a new DataFrame for plotting
data['subjChr'] = data['CHROM'].apply(get_chrom_number)
data.dropna(subset=['subjChr'], inplace=True)

# Read the centromere file
chrcen_data = []
with open(centromere_file) as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

# Read the genes file
genes_data = []
with open(genes_file) as f:
    for line in f:
        chrom, pos, gene = line.split()
        genes_data.append((chrom, int(pos), gene))

# Read the chromosome end file
end_data = []
with open(end_file) as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom, int(pos)))

# Define marker shapes and colors for different donors
marker_shapes = {
    'D1T4_calling.vcf': 'o',  # Circle
    'D2T4_calling.vcf': 'o',  # Square
    'D4T4_calling.vcf': 'o',  # Triangle
    'D6T4b_calling.vcf': '^',  # Right Triangle
    'D6T4a_merged.vcf': '^'   # Left Triangle
}

marker_colors = {
    'D1T4_calling.vcf': '#d62728',  # Red
    'D2T4_calling.vcf': '#9467bd',  # Purple
    'D4T4_calling.vcf': '#8c564b',  # Brown
    'D6T4b_calling.vcf': '#e377c2',  # Pink
    'D6T4a_merged.vcf': '#7f7f7f'   # Gray
}

# Define custom labels for the legend
custom_labels = {
    'D1T4_calling.vcf': 'Donor 1',
    'D2T4_calling.vcf': 'Donor 2',
    'D4T4_calling.vcf': 'Donor 4',
    'D6T4b_calling.vcf': 'Donor 6b',
    'D6T4a_merged.vcf': 'Donor 6a'
}

# Define the plot
fig, ax = plt.subplots(figsize=(12, 8))
plt.tight_layout(pad=4)

# Plot each sample with a different color and marker
x_offset_multiplier = 2  # Adjust this multiplier to control space between chromosomes
for sample in data['Sample'].unique():
    sample_data = data[data['Sample'] == sample]

    # Get the marker shape and color for this sample
    marker = marker_shapes.get(sample, 'o')  # Default to circle if sample not in dictionary
    color = marker_colors.get(sample, '#000000')  # Default to black if sample not in dictionary

    # Plot with the label for the legend
    ax.scatter(sample_data['subjChr'] * x_offset_multiplier, sample_data['POS'],
               label=sample, alpha=0.5, marker=marker, color=color)

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data):
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

        # Add ellipse for the first half of the chromosome
        ax.add_patch(Ellipse((chrom_num * x_offset_multiplier, centromere_position / 2),
                             0.35, centromere_position, edgecolor='black', linewidth=0.4, fill=False, zorder=20))

        # Add ellipse for the second half of the chromosome
        ax.add_patch(Ellipse((chrom_num * x_offset_multiplier, (end + centromere_position) / 2),
                             0.35, end - centromere_position, edgecolor='black', linewidth=0.4, fill=False, zorder=20))

# Add centromere diamonds
for chrom, cen in chrcen_data:
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        ax.plot([chrom_num * x_offset_multiplier], [cen], marker='D', markersize=3.75,
                markerfacecolor="grey", color='grey', zorder=21)

# Add the genes data as text labels with lines connecting to their chromosome positions
for chrom, pos, gene in genes_data:
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        x_pos = chrom_num * x_offset_multiplier
        y_pos = pos
        ax.annotate(gene, xy=(x_pos, y_pos), xytext=(x_pos + 0.25, y_pos),
                    textcoords='data', ha='left', va='center', fontsize=11, rotation=90,
                    fontweight='bold', arrowprops=dict(arrowstyle="->", color='black', lw=0.8))

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Retroseq Human Swim Up T4', fontsize=20, fontweight='bold')
chrom_labels = [i * x_offset_multiplier for i in range(1, 25)]
ax.set_xticks(chrom_labels)
ax.set_xticklabels([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'], fontsize=8.5, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend with custom labels in the specified order
legend_elements = [
    Line2D([0], [0], marker=marker_shapes.get(sample, 'o'), color='w', markersize=10,
           markerfacecolor=marker_colors.get(sample, '#000000'), label=custom_labels.get(sample, sample))
    for sample in custom_labels.keys()  # Use keys of custom_labels for order
]

plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'DMAG_SU_T4_retroseq.png'  # Desired filename for the output plot
plt.savefig(output_file, dpi=600)
plt.show()
