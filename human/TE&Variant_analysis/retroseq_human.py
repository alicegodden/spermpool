#Title: Retroseq plotting
#Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import numpy as np

def extract_chrom_pos_from_vcf(vcf_file):
    # Initialize dictionaries for chromosome positions and TE categories
    chrom_pos_dict = {str(i): [] for i in range(1, 23)}  # Map '1'-'22' to 1-22
    chrom_pos_dict['X'] = []  # Add X
    chrom_pos_dict['Y'] = []  # Add Y
    te_category_dict = {str(i): [] for i in range(1, 23)}
    te_category_dict['X'] = []  # Add X
    te_category_dict['Y'] = []  # Add Y

    # Read the VCF file using pysam
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            # Replace 'chr' if present in chromosome name for consistency
            chromosome = record.chrom.replace("chr", "")
            if chromosome in chrom_pos_dict:
                position = record.pos
                info_field = record.info.get('MEINFO', '')

                # Ensure info_field is a string for processing
                if isinstance(info_field, tuple):
                    info_field = info_field[0]
                elif not isinstance(info_field, str):
                    info_field = str(info_field)

                # Determine the type of transposable element
                if info_field.strip():
                    te_type = info_field.split()[0]
                    te_category = (
                        'DNA' if 'DNA' in te_type else
                        'LINE' if 'LINE' in te_type else
                        'LTR' if 'LTR' in te_type else
                        'SINE' if 'SINE' in te_type else
                        'RC' if 'RC' in te_type else
                        'Satellite' if 'SATELLITE' in te_type else
                        'Unknown'
                    )
                else:
                    te_category = 'Unknown'

                # Append the position and TE category to the dictionaries
                chrom_pos_dict[chromosome].append(position)
                te_category_dict[chromosome].append(te_category)
    return chrom_pos_dict, te_category_dict


# Replace with your VCF file path
vcf_file = 'hapselDK_raw_merged_chrom_win100_gq35_fl8.vcf'
chrom_pos_dict, te_category_dict = extract_chrom_pos_from_vcf(vcf_file)

# Read the centromere data
chrcen_data = []
with open('centromere_grch38.txt') as f:  # Ensure correct file extension
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom.replace("chr", ""), int(pos)))  # Strip 'chr'

# Read the genes data
genes_data = []
with open('GRCh38_genes.txt') as f:  # Ensure correct file extension
    for line in f:
        chrom, pos, gene = line.split()
        genes_data.append((chrom.replace("chr", ""), int(pos), gene))

# Read chromosome lengths
end_data = []
with open('GRCh38_chrom_length.txt') as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom.replace("chr", ""), int(pos)))

# Color palette for TE types
flare_palette = sns.color_palette("twilight")
te_colors = {
    'DNA': flare_palette[0],
    'LINE': flare_palette[2],
    'LTR': flare_palette[3],
    'RC': flare_palette[4],
    'Satellite': flare_palette[5],
    'SINE': flare_palette[1],
    'Unknown': 'grey'
}

# Create a mapping for chromosomes to numeric values
chromosome_map = {str(i): i for i in range(1, 23)}  # Map '1'-'22' to 1-22
chromosome_map['X'] = 23  # Assign X a numeric value (23)
chromosome_map['Y'] = 24  # Assign Y a numeric value (24)

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Plotting TE insertion points with alpha transparency
for chromosome, positions in chrom_pos_dict.items():
    # Check if chromosome is in the mapping
    if chromosome in chromosome_map:
        numeric_chromosome = chromosome_map[chromosome]  # Use mapping for numeric chromosome
        categories = te_category_dict[chromosome]
        colors = [te_colors[category] for category in categories]
        ax.scatter(
            [numeric_chromosome] * len(positions),
            positions,
            color=colors,
            alpha=0.5,
            label='Chromosome ' + chromosome
        )
    else:
        print(f"Warning: Chromosome '{chromosome}' not found in mapping. Skipping.")

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data):
    # Use the chromosome mapping to get the numeric value
    if chrom in chromosome_map:
        numeric_chrom = chromosome_map[chrom]
        centromere_position = next((cen for chr_name, cen in chrcen_data if chr_name == chrom), None)

        if centromere_position is not None:
            # First half of the chromosome (from 0 to centromere)
            ax.add_patch(
                Ellipse(
                    (numeric_chrom, centromere_position / 2),
                    0.27,
                    centromere_position,
                    edgecolor='black',
                    linewidth=0.4,
                    fill=False,
                    zorder=20
                )
            )

            # Second half of the chromosome (from centromere to end)
            ax.add_patch(
                Ellipse(
                    (numeric_chrom, (end + centromere_position) / 2),
                    0.27,
                    end - centromere_position,
                    edgecolor='black',
                    linewidth=0.4,
                    fill=False,
                    zorder=20
                )
            )
    else:
        print(f"Warning: Chromosome '{chrom}' not found in mapping for ellipses. Skipping.")

# Add centromere markers
for chrom, cen in chrcen_data:
    ax.plot([chromosome_map[chrom]], [cen], marker='D', markersize=3.75, markerfacecolor="grey", color='grey',
            zorder=21)

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Non-reference TE Insertion mutations Unselected', fontsize=20, fontweight='bold')

# Create a list of chromosome labels with 'chr' prefix
chromosome_labels = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
plt.xticks(range(1, 25), chromosome_labels, fontsize=8, fontweight='bold')  # Updated to use chromosome_labels
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend for both plots
legend_elements = [
                      Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, alpha=0.5,
                             label=te_type)
                      for te_type, color in te_colors.items() if te_type != 'Unknown'
                  ] + [
                      Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
                  ]
plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'Haploid_retroseq_DMAG_gq30_fl8_unselRAW.png'  # Replace with your desired filename
plt.savefig(output_file, dpi=600)
plt.show()
