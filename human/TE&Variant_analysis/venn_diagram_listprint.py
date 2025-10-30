# Title: Venn diagram
# Author: Dr. Alice M. Godden
import matplotlib.pyplot as plt
from venn import venn


# Function to read gene lists from text files
def read_gene_list(file_path):
    with open(file_path, 'r') as file:
        genes = set(line.strip() for line in file)
    return genes


# File paths for your input files
files = {
    'Zebrafish_Sperm_human_orthologues': 'fish2hum_sort',  # Replace with actual file paths
    'T2_MethylCellulose_Exp1': 'exp1_MC_donors1_4',  # Replace with actual file paths
    'T4_SwimUp_Exp2': 'exp2_swimup_donors_5_9_t4',  # Replace with actual file paths
    #'T8_SwimUp_Exp2': 'T0vT8',
#'T12_SwimUp_Exp2': 'T0vT12',
#'T24_SwimUp_Exp2': 'T0vT24',
#'T48_SwimUp_Exp2': 'T0vT48',
}

# Read gene lists into a dictionary
gene_sets = {name: read_gene_list(path) for name, path in files.items()}


# Function to calculate unique genes for each group
def get_unique_genes_for_all_groups(gene_sets):
    unique_genes = {}

    # For each group, subtract the union of all other groups
    for group_name in gene_sets:
        others = set.union(*(gene_sets[name] for name in gene_sets if name != group_name))
        unique_genes[group_name] = gene_sets[group_name] - others

    return unique_genes


# Get unique genes for each group
unique_genes = get_unique_genes_for_all_groups(gene_sets)

# Find genes shared by all groups
shared_genes = set.intersection(*gene_sets.values())

# Save all unique genes and shared genes into a single file with headers
with open('unique_genes_all_Fish_T2_T4.txt', 'w') as file:
    # Write the shared genes
    file.write("Genes shared by all groups:\n")
    for gene in sorted(shared_genes):
        file.write(f"{gene}\n")

    file.write("\n\n")

    # Write the unique genes for each group
    for group, unique_genes_list in unique_genes.items():
        # Write the header
        file.write(f"Unique genes in '{group}' that are not in any other group:\n")

        # Write the genes
        for gene in sorted(unique_genes_list):
            file.write(f"{gene}\n")

        # Add a separator between different group sections
        file.write("\n\n")

        print(f"\nSaved the list of unique genes for '{group}' to 'unique_genes_all_groups.txt'.\n")

# Create a dictionary in the format required by the venn library
venn_data = {key: value for key, value in gene_sets.items()}

# Plot the Venn diagram using the venn library
plt.figure(figsize=(14, 14))
venn_diagram = venn(venn_data)

# Customize for publication quality
plt.title('Zebrafish vs MethylCellulose T2 Exp 1 vs Swim-Up T4 Exp 2', fontsize=14, fontweight='bold')

# Set text size and bold for all text in the Venn diagram
for text in plt.gca().texts:  # Access all text objects in the current Axes
    text.set_fontweight('bold')
    text.set_fontsize(20)  # Change the size to your desired font size

# Make legend text bold
legend = plt.gca().get_legend()
if legend:
    for text in legend.get_texts():
        text.set_fontweight('bold')
        text.set_fontsize(16)  # Change the size to your desired font size

plt.tight_layout()
# Save the figure in a publication-quality format
plt.savefig('Venn_Diagram_Genes_FishvsExp1vExp2.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# Print the venn_data dictionary
print("\nVenn data (all genes in each set):")
print(venn_data)
