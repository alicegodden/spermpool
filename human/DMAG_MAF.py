# Title: Plotting MAF on chromosomal plots
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import numpy as np

# Define file paths
input_file = 'all_filt_males_MafO_0.95.csv'  # Replace with your CSV file path

# Read the input CSV file
data = pd.read_csv(input_file)

# Check for NaN values in the 'Sample' column and handle them
if data['Sample'].isnull().any():
    print("Found NaN values in 'sample' column. Removing rows with NaN values.")
    data = data.dropna(subset=['Sample'])

# Assign colors to samples
samples = data['Sample'].unique()
palette = sns.color_palette("Blues", len(samples))
sample_colors = {sample: color for sample, color in zip(samples, palette)}

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Read the chrcen.txt file
chrcen_data = []
with open('chrcen.txt') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

# Read the genes.txt file
genes_data = []
with open('zfish_DMAG_genes') as f:
    for line in f:
        chrom, pos, gene = line.split()
        genes_data.append((chrom, int(pos), gene))

# Read the chrom_end.txt file
end_data = []
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom, int(pos)))

# Add scatter plot for each sample
for sample in samples:
    sample_data = data[data['Sample'] == sample]
    if 'subjChr' in sample_data.columns and 'subjStart' in sample_data.columns:
        ax.scatter(sample_data['subjChr'].astype(int), sample_data['subjStart'], color=sample_colors[sample],
                   label=sample, alpha=0.5)

# Add ellipses for chromosomes
for chrom, end in end_data:
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.4,
                fill=False, zorder=20))

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.4, fill=False, zorder=20))


# Define the function to map chromosome names to numeric values
def get_chrom_number(chrom):
    try:
        return int(chrom)
    except ValueError:
        # Handle non-numeric chromosome names if needed
        chrom_map = {'X': 25, 'Y': 26}  # Example for non-numeric chromosomes
        return chrom_map.get(chrom, None)


# Function to prevent overlapping labels by adjusting their positions
def adjust_label_position(y_pos, used_positions, min_distance=50):
    """ Adjusts the y-position of a label to avoid overlapping with other labels """
    for used_y in used_positions:
        if abs(y_pos - used_y) < min_distance:
            y_pos += min_distance  # Shift the label down if too close to an existing label
    used_positions.append(y_pos)
    return y_pos


# Track used positions to avoid label overlap
used_positions = []

# Add gene annotations with arrows and rotation of 90 degrees, place labels on the left
x_offset_multiplier = 0.27  # Use as a factor to offset the gene label from the chromosome

for chrom, pos, gene in genes_data:
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        # Offset for labels to the left of the chromosome ellipse
        x_offset = 0.25  # Adjust this value as needed
        x_pos = chrom_num  # Chromosome position
        y_pos = pos

        # Adjust y-position to avoid overlap with other labels
        y_pos = adjust_label_position(y_pos, used_positions, min_distance=50)  # Set min distance between labels

        # Annotate with an offset to the left of the chromosome position
        ax.annotate(
            gene,
            xy=(x_pos, y_pos),  # Position of the gene
            xytext=(x_pos - x_offset, y_pos),  # Position of the text with offset to the left
            textcoords='data',
            ha='right',  # Align text to the right (to move it to the left of the point)
            va='center',
            fontsize=6.75,
            rotation=90,  # Rotate labels 90 degrees
            fontweight='bold',
            arrowprops=dict(
                arrowstyle="->",  # Adds the line joining the label and the position
                color='black',
                lw=0.8,
            )
        )

# Add centromere diamonds
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=3.75, markerfacecolor="grey", color='grey', zorder=21)


# Function to read dopes data (chromatin accessibility regions)
def read_dopes_file(dopes_file):
    dopes_data = []
    with open(dopes_file) as f:
        next(f)  # Skip header line
        for line in f:
            chrom, chromStart, chromEnd, name, score, strand = line.split()
            chrom = chrom.replace("chr", "")
            dopes_data.append((chrom, int(chromStart), int(chromEnd), name, int(score), strand))
    return dopes_data


# Add chromatin accessibility data (open regions in blue) next to the ellipse on the right side
for chrom, chromStart, chromEnd, name, score, strand in read_dopes_file('daniocode_hub_280355_dopes_all.txt'):
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        ax.hlines(y=[chromStart, chromEnd], xmin=chrom_num + 0.15, xmax=chrom_num + 0.25, color='darkblue',
                  linewidth=0.5, zorder=30)

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Unselected Zebrafish Maf_O > 0.95 cov > 200', fontsize=20, fontweight='bold')
plt.xticks(range(1, 26), fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=sample)
    for sample, color in sample_colors.items()
] + [
    Line2D([0], [0], color='darkblue', lw=2, label='Open regions'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
]
plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'allmale_zfish_Maf0.95_cov200_2409.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)
plt.show()
