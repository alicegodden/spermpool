import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

# Define file paths
input_file = 'Males_zebrafish.csv'  # Replace with your CSV file path
# Read the input CSV file
data = pd.read_csv(input_file)

# Assign colors to samples
samples = data['Sample'].unique()
palette = sns.color_palette("Blues", len(samples))
sample_colors = {sample: color for sample, color in zip(samples, palette)}

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Plot each sample with a different color

# Add ellipses for chromosomes based on length
# Read the chrcen.txt file
chrcen_data = []
with open('chrcen.txt') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))
# Read the genes.txt file
genes_data = []
with open('genes.txt') as f:
    for line in f:
        chrom, pos, gene = line.split()
        genes_data.append((chrom, int(pos), gene))
# Read the chrom_end.txt file
end_data = []
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom, int(pos)))


# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))



# Add the genes data as text labels with marks
for sample in samples:
    sample_data = data[data['Sample'] == sample]
    ax.scatter(sample_data['subjChr'].astype(int), sample_data['subjStart'], color=sample_colors[sample], label=sample,
               alpha=0.6)

for chrom, pos, gene in genes_data:
    x, y = int(chrom), pos
    ax.annotate(gene, xy=(x, y), ha='center', va='center', fontsize=8, fontweight='bold',
                arrowprops=dict(arrowstyle="<-", color='red', facecolor='none', linewidth=0))

# Add ellipses for each half of the chromosome
ellipse_handles = []
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.4,
                fill=False, zorder=20))

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.4, fill=False, zorder=20))

# Add centromere diamonds
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=3.75, markerfacecolor="grey", color='grey', zorder=21)

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Gene Annotations by Sample', fontsize=20, fontweight='bold')
plt.xticks(range(1, 26), fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend
legend_elements = [
                      Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=sample)
                      for sample, color in sample_colors.items()
                  ] + [
                      Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')]
plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'gene_annotations_by_sample.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)
plt.show()
