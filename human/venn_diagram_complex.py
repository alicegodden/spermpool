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
    'Fish': 'fish2hum_sort',
    'T2': 'T0vT2',
    'T4': 'T0vT4',
    'T48': 'T0vT48'
}

# Read gene lists into a dictionary
gene_sets = {name: read_gene_list(path) for name, path in files.items()}

# Find the overlapping items shared across all lists
shared_items = set.intersection(*gene_sets.values())

# Write the shared items to a text file
with open('shared.txt', 'w') as file:
    for item in sorted(shared_items):
        file.write(f"{item}\n")

# Create a dictionary in the format required by the venn library
venn_data = {key: value for key, value in gene_sets.items()}

# Plot the Venn diagram using the venn library
plt.figure(figsize=(10, 8))
venn_diagram = venn(venn_data)

# Customize for publication quality
plt.title('Venn Diagram of fishy human spermy genes', fontsize=14, fontweight='bold')

# Set text size and bold for all text in the Venn diagram
for text in plt.gca().texts:  # Access all text objects in the current Axes
    text.set_fontweight('bold')
    text.set_fontsize(20)  # Change the size to your desired font size (e.g., 12)

# Make legend text bold
legend = plt.gca().get_legend()
if legend:
    for text in legend.get_texts():
        text.set_fontweight('bold')
        text.set_fontsize(16)  # Change the size to your desired font size (e.g., 12)

plt.tight_layout()

# Save the figure in a publication-quality format
plt.savefig('Venn_Diagram.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
