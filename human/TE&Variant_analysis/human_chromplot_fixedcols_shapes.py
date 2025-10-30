# Title: Chromosomal plots of variants
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re  # Import the re module for regular expressions
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

# Define file paths
input_file = 'T0xT2_cov200.csv'  # Replace with your CSV file path

# Read the input CSV file
data = pd.read_csv(input_file)

# Check for NaN values in the 'Sample' column and handle them
if data['Sample'].isnull().any():
    print("Found NaN values in 'Sample' column. Removing rows with NaN values.")
    data = data.dropna(subset=['Sample'])

# Define a specific color palette for the donors
donor_palette = {
    'Donor1': '#1f77b4',  # Blue
    'Donor2': '#ff7f0e',  # Orange
    'Donor3': '#2ca02c',  # Green
    'Donor4': '#d62728',  # Red
    'Donor5': '#9467bd',  # Purple
    'Donor6': '#8c564b',  # Brown
    #'Donor7': '#e377c2',  # Pink
    #'Donor8': '#7f7f7f',  # Gray
    #'Donor9': '#bcbd22'  # Olive
}

# Get the colors for each sample
sample_colors = {sample: donor_palette.get(sample, '#000000') for sample in data['Sample'].unique()}

# Read the centromere file
chrcen_data = []
with open('centromere_grch38') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

# Read the genes file
genes_data = []
with open('GRCh38_genes') as f:
    for line in f:
        chrom, pos, gene = line.split()
        genes_data.append((chrom, int(pos), gene))

# Read the chromosome end file
end_data = []
with open('GRCh38_chrom_length.txt') as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom, int(pos)))

# Function to map chromosomes
def get_chrom_number(chrom):
    if chrom == "chrX":
        return 23
    elif chrom == "chrY":
        return 24
    else:
        try:
            chrom_num = int(chrom.replace("chr", ""))
            if 1 <= chrom_num <= 22:
                return chrom_num
        except ValueError:
            pass
    return None  # Exclude any other chromosomes

# Define the list of chromosomes in the correct order
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# Filter data to include only chromosomes 1-22, X, and Y
data['subjChr'] = data['subjChr'].apply(get_chrom_number)
data = data.dropna(subset=['subjChr'])

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))
plt.tight_layout(pad=4)

# Plot each sample with a different color and marker
x_offset_multiplier = 2  # Adjust this multiplier to control space between chromosomes
for sample in data['Sample'].unique():
    sample_data = data[data['Sample'] == sample]

    # Determine the marker based on the donor number
    donor_number = int(re.search(r'\d+', sample).group())  # Extract the donor number from the sample name
    if 5 <= donor_number <= 9:
        marker = '^'  # Use '^' marker (triangle) for donors 5-9
        marker_size = 80
    else:
        marker = 'o'  # Use 'o' marker (circle) for donors 1-4
        marker_size = 60

    # Plot the data with the corresponding color and marker
    ax.scatter(sample_data['subjChr'] * x_offset_multiplier, sample_data['subjStart'],
               color=sample_colors[sample], marker=marker, label=sample, alpha=0.95)

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data):
    chrom_num = get_chrom_number(chrom)
    if chrom_num:
        centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

        # Add ellipse for the first half of the chromosome (from 0 to centromere)
        ax.add_patch(
            Ellipse((chrom_num * x_offset_multiplier, centromere_position / 2), 0.35, centromere_position,
                    edgecolor='black', linewidth=0.4, fill=False, zorder=20))

        # Add ellipse for the second half of the chromosome (from centromere to end)
        ax.add_patch(
            Ellipse((chrom_num * x_offset_multiplier, (end + centromere_position) / 2), 0.35, end - centromere_position,
                    edgecolor='black', linewidth=0.4, fill=False, zorder=20))

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
        # Offset for labels to the right of the chromosome ellipse
        x_offset = 0.25  # Adjust this value as needed
        x_pos = chrom_num * x_offset_multiplier  # Use chrom_num directly from the mapping
        y_pos = pos

        # Annotate with an offset to the right of the chromosome position
        ax.annotate(
            gene,
            xy=(x_pos, y_pos),  # Position of the gene
            xytext=(x_pos + x_offset * 2, y_pos),  # Position of the text with offset
            textcoords='data',
            ha='left',  # Align text to the left (to move it to the right of the point)
            va='center',
            fontsize=11,
            rotation=90,
            fontweight='bold',
            arrowprops=dict(
                arrowstyle="->",  # Adds the line joining the label and the position
                color='black',
                lw=0.8,
            )
        )

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Human Experiment 1 (Methyl Cellulose)', fontsize=20, fontweight='bold')
chrom_labels = [i * x_offset_multiplier for i in range(1, 25)]
ax.set_xticks(chrom_labels)
ax.set_xticklabels([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'], fontsize=8.5, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=sample_colors[sample], markersize=10, label=sample)
    for sample, color in sample_colors.items() if int(re.search(r'\d+', sample).group()) <= 4
] + [
    Line2D([0], [0], marker='^', color='w', markerfacecolor=sample_colors[sample], markersize=10, label=sample)
    for sample, color in sample_colors.items() if 5 <= int(re.search(r'\d+', sample).group()) <= 9
] + [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
]

plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot

output_file = 'DMAG_T0vt2_EXP1_title.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)
plt.show()
