Title: TE counts per chromosome
Subtitle: Vcf file retroseq output as input
Author: Dr. Alice M. Godden



import matplotlib.pyplot as plt
import pysam
import seaborn as sns

def extract_chrom_pos_from_vcf(vcf_file):
    chrom_pos_dict = {str(i): [] for i in range(1, 26)}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            if chromosome in chrom_pos_dict:
                position = record.pos
                chrom_pos_dict[chromosome].append(position)
    return chrom_pos_dict

# Replace with your VCF file paths
vcf_files = [
    'Non_minusNon_nmdups_win100_chrom_gq1000_fl8.vcf',
    'Finclip_minusNon_nmdups_win100_chrom_gq1000_fl8.vcf',
    'Central_minusNon_nmdups_win100_chrom_gq1000_fl8.vcf',
    'Outer_minusNon_nmdups_win100_chrom_gq1000_fl8.vcf'
]
# Corresponding sample names
sample_names = [
    'Non-selected',
    'Fin-clip',
    'Central',
    'Outer'
]

# Store counts for each VCF file
all_te_counts = []

# Process each VCF file
for vcf_file in vcf_files:
    chrom_pos_dict = extract_chrom_pos_from_vcf(vcf_file)
    te_counts = {chrom: len(positions) for chrom, positions in chrom_pos_dict.items()}
    all_te_counts.append(te_counts)

# Prepare data for plotting
chromosomes = list(map(str, range(1, 26)))
data = {sample_name: [te_counts.get(chrom, 0) for chrom in chromosomes] for sample_name, te_counts in zip(sample_names, all_te_counts)}

# Color palette for the lines
rocket_palette = sns.color_palette("mako", len(sample_names))

# Plotting
fig, ax = plt.subplots(figsize=(15, 8))

for sample_name, color in zip(sample_names, rocket_palette):
    ax.plot(chromosomes, data[sample_name], label=sample_name, color=color, marker='o')

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Counts of TEs', fontsize=18, fontweight='bold')
plt.title('TEs per Chromosome from Retroseq', fontsize=20, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.legend(title='Samples', fontsize=12, title_fontsize=14, loc='upper right')

# Save and display the plot
output_file = 'TE_counts_per_chromosome_spermpool.png'
plt.savefig(output_file, dpi=600)
plt.show()
