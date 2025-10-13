# Title: Heatmap plotting z-score normalised 99th percentile data for polygenic risk scores
# Author: Dr. Alice M. Godden

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- SET INTERACTIVE PLOTTING MODE (optional) ---
plt.ion()

# --- HELPER FUNCTION TO CLEAN TRAIT LABELS ---
def clean_label(label):
    """
    Cleans the filename to extract a more readable trait label.
    """
    label_str = str(label)
    return label_str.replace("_filtered.tsv.bgz", "").split("-both_sexes")[0]

# --- MAIN DATA PROCESSING FUNCTION FOR AN INDIVIDUAL FILE ---
def process_individual_file(filepath, phenocodes_df):
    """
    Reads an individual's score difference file, cleans labels, and merges phenocodes.
    ASSUMPTION: File contains match_file and score_difference for one individual.
    """
    # Extract Individual ID from the filename provided in the input_files list
    # e.g., "PGSscores_difference_D1T4-D1T0.txt" -> "D1T4-D1T0"
    individual_id = os.path.basename(filepath).replace("PGSscores_difference_", "").replace(".txt", "")

    try:
        df_individual = pd.read_csv(filepath, sep='\t')
        if df_individual.empty:
            print(f"  - Warning: Individual file {filepath} is empty. Skipping.")
            return None

        required_cols = ['match_file', 'score_difference']
        if not all(col in df_individual.columns for col in required_cols):
            print(f"  - Error: Missing required columns ({required_cols}) in individual file {filepath}. Skipping.")
            print(f"    - Found columns: {df_individual.columns.tolist()}")
            return None

        # Clean trait labels
        df_individual["label"] = df_individual["match_file"].apply(clean_label)
        df_individual["filename_clean"] = df_individual["match_file"].str.replace("_filtered.tsv.bgz", ".tsv.bgz")

        # Merge with phenocodes for descriptive labels
        merged_df = df_individual.merge(phenocodes_df, left_on="filename_clean", right_on="filename", how="left")
        merged_df["label"] = merged_df.apply(
            lambda row: f"{row['label']} → {row['description']}" if pd.notna(row["description"]) else row["label"],
            axis=1
        )

        merged_df['score_difference'] = pd.to_numeric(merged_df['score_difference'], errors='coerce')
        merged_df_cleaned = merged_df.dropna(subset=['score_difference', 'label'])

        if merged_df_cleaned.empty:
            print(f"  - Warning: No valid score_difference data in '{filepath}' after cleaning. Skipping.")
            return None

        # --- Handle Duplicate Trait Labels for a single individual ---
        if merged_df_cleaned['label'].duplicated().any():
            print(f"    - Warning: Duplicate trait labels found in {filepath}. Taking mean score_difference for duplicates.")
            merged_df_cleaned = merged_df_cleaned.groupby('label')['score_difference'].mean().reset_index()

        # Ensure 'label' column is unique before setting as index
        if not merged_df_cleaned['label'].is_unique:
             print(f"    - Critical Error: Trait labels are still not unique in {filepath} after mean aggregation. This should not happen. Skipping.")
             return None

        return merged_df_cleaned[['label', 'score_difference']].set_index('label').rename(
            columns={'score_difference': individual_id} # Use the extracted individual_id as column name
        )

    except Exception as e:
        print(f"  - An error occurred while reading individual file {filepath}: {e}")
        return None

# --- CONFIGURATION ---
# List of individual files you want to include in the heatmap

# SU (commented out)
input_files = [
#    "PGSscores_difference_D1T4-D1T0.txt",
    "PGSscores_difference_D2T4-D2T0.txt",
#    "PGSscores_difference_D4T4-D4T0.txt",
#    "PGSscores_difference_D6aT4-D6T0a.txt",
#    "PGSscores_difference_D6bT4-D6bT0.txt",
#]

# MC (active)
#input_files = [
#    "PGSscores_difference_MC_M8_C-O.txt",
#    "PGSscores_difference_MC_M11_C-O.txt",
#   "PGSscores_difference_MC_M12_C-O.txt",
#    "PGSscores_difference_MC_M13_C-O.txt",
]
# Ensure the script can find these files. If they are not in the same directory
# as the script, you might need to provide their full paths or a relative path.
# Example: input_files = ["/path/to/my/data/PGSscores_difference_D1T4-D1T0.txt", ...]

phenocodes_path = "phenocodes" # Ensure this path is correct relative to where you run the script, or an absolute path
output_heatmap_filename = "PGSscore_heatmap_Individuals_Normalized_AnnotatedRaw_20_SU_D1.png" # Updated output filename
PERCENTILE_THRESHOLD = 99 # Filter by top N percentile of overall absolute differences across all individuals/traits
MAX_TRAITS_TO_PLOT = 20 # Keep a cap for readability

# --- LOAD PHENOCODES ---
print(f"\n--- Loading phenocodes from: {phenocodes_path} ---")
try:
    phenos = pd.read_csv(phenocodes_path, sep='\t')
    print(f"  - Successfully loaded phenocodes. Shape: {phenos.shape}")
    if 'filename' not in phenos.columns or 'description' not in phenos.columns:
        print(f"  - Error: 'filename' or 'description' column missing in '{phenocodes_path}'. Please check phenocodes file.")
        exit(1)
except FileNotFoundError:
    print(f"Error: Phenocodes file not found at '{phenocodes_path}'. Please check the path.")
    exit(1)
except Exception as e:
    print(f"Error loading phenocodes: {e}")
    exit(1)

# --- PROCESS ALL INDIVIDUAL FILES ---
all_individuals_dfs = []
print(f"\n--- Processing specified individual files ---")

if not input_files:
    print(f"Error: No input files specified in the 'input_files' list. Exiting.")
    exit(1)

for filepath in input_files:
    if os.path.exists(filepath):
        print(f"  - Processing '{filepath}'...")
        processed_df = process_individual_file(filepath, phenos)
        if processed_df is not None:
            all_individuals_dfs.append(processed_df)
    else:
        print(f"Warning: Input file not found: '{filepath}'. Skipping.")


if not all_individuals_dfs:
    print("No valid individual data found from input files. Exiting.")
    exit(0)

# Combine all individual DataFrames into a single heatmap-ready DataFrame
# This will result in traits as rows, individual IDs as columns
heatmap_data_raw = pd.concat(all_individuals_dfs, axis=1)

print(f"\n--- Combined raw data from all individuals for heatmap ---")
print(f"  - Final heatmap_data_raw DataFrame shape (before filtering): {heatmap_data_raw.shape}")
print(f"  - Head of heatmap_data_raw (before filtering):")
print(heatmap_data_raw.head())
print(f"  - Number of NaN values in heatmap_data_raw (before filtering): {heatmap_data_raw.isnull().sum().sum()}")


# --- TRAIT SELECTION (Percentile Filtering based on overall max absolute difference) ---
print(f"\n--- Trait Selection (Percentile Filtering) ---")

# Calculate the maximum absolute score difference for each trait across all individuals
# Use .replace([np.inf, -np.inf], np.nan) to handle any infinities first, then .max()
max_abs_diff_per_trait = heatmap_data_raw.replace([np.inf, -np.inf], np.nan).stack().abs().groupby(level=0).max().sort_values(ascending=False)

# Remove any NaNs from the max_abs_diff_per_trait series itself, as they can mess with quantile calculation
max_abs_diff_per_trait = max_abs_diff_per_trait.dropna()


print(f"  - Max absolute difference per trait (top 5):")
print(max_abs_diff_per_trait.head())

total_unique_traits = len(max_abs_diff_per_trait)

# Calculate the threshold value for the PERCENTILE_THRESHOLD
if total_unique_traits > 0:
    # Ensure quantile calculation is robust even if most values are 0.0 or very small
    # If percentile_value ends up being 0.0 due to many 0s, it's fine, it means all non-zero will qualify.
    percentile_value = max_abs_diff_per_trait.quantile(1 - (PERCENTILE_THRESHOLD / 100.0))

    # Traits where max absolute difference is strictly greater than 0, if percentile_value is 0.0
    if percentile_value == 0.0:
        qualifying_traits_series = max_abs_diff_per_trait[max_abs_diff_per_trait > 0]
    else:
        qualifying_traits_series = max_abs_diff_per_trait[max_abs_diff_per_trait >= percentile_value]

    selected_labels = qualifying_traits_series.head(MAX_TRAITS_TO_PLOT).index.tolist()

    print(f"  - Total unique traits available: {total_unique_traits}")
    print(f"  - Threshold value for {PERCENTILE_THRESHOLD}th percentile (absolute difference): {percentile_value:.4f}")
    print(f"  - Number of traits meeting or exceeding the {PERCENTILE_THRESHOLD}th percentile: {len(qualifying_traits_series)}")

    if not selected_labels and total_unique_traits > 0:
        print(f"    (No traits strictly above {percentile_value:.4f} or MAX_TRAITS_TO_PLOT is 0. Selecting top 1 trait as fallback.)")
        selected_labels = max_abs_diff_per_trait.head(1).index.tolist()
else:
    selected_labels = []

print(f"  - Final number of labels selected for plotting (max {MAX_TRAITS_TO_PLOT}): {len(selected_labels)}")

if not selected_labels:
    print("No traits met the selection criteria for plotting. Cannot create heatmap. Exiting.")
    exit(0)

# Filter the raw data to include only these selected labels
# This filtered raw data will be used for annotations
filtered_heatmap_data_raw_for_annot = heatmap_data_raw.loc[selected_labels].copy()


# --- Z-SCORE NORMALIZATION FOR PLOTTING (Within Each Individual) ---
print(f"\n--- Normalizing score differences within each individual for heatmap ---")
normalized_heatmap_data = filtered_heatmap_data_raw_for_annot.copy()

# MODIFIED: More robust handling of zero standard deviation:
# If std is 0, all values are the same; Z-score should be 0.
# If a column is all NaN, it will remain NaN.
def zscore_normalize(series):
    if series.std() == 0:
        return series * 0.0 # Return a series of zeros if std is zero
    else:
        return (series - series.mean()) / series.std()

normalized_heatmap_data = normalized_heatmap_data.apply(zscore_normalize, axis=0)

# IMPORTANT: Ensure that the indices (traits) and columns (individuals)
# of normalized_heatmap_data and filtered_heatmap_data_raw_for_annot
# remain aligned after any drops due to NaNs.
# We'll apply the same dropna operations to both.

# Drop NaNs from normalized data (this will determine the final shape for plotting)
# Only drop columns that are ALL NaN (individuals with no valid data after normalization/selection)
normalized_heatmap_data.dropna(axis=1, how='all', inplace=True)
# Only drop rows that are ALL NaN (traits with no valid data after normalization/selection for any remaining individual)
normalized_heatmap_data.dropna(axis=0, how='all', inplace=True)


# Now, align the raw data for annotation to the same shape as normalized data
# Only keep columns and rows that are still present in normalized_heatmap_data
# This ensures annot and data have the same shape.
filtered_heatmap_data_raw_for_annot = filtered_heatmap_data_raw_for_annot.loc[
    normalized_heatmap_data.index, normalized_heatmap_data.columns
]


if normalized_heatmap_data.empty:
    print("No data remaining after normalization. Cannot create heatmap. Exiting.")
    exit(0)

# Sort traits by their overall normalized score for consistent plotting order
# This ensures traits with generally higher (or lower) normalized scores appear together
normalized_heatmap_data['overall_normalized_score'] = normalized_heatmap_data.mean(axis=1)
normalized_heatmap_data = normalized_heatmap_data.sort_values(
    by='overall_normalized_score', ascending=False).drop(columns='overall_normalized_score')

# Re-align filtered_heatmap_data_raw_for_annot to the final sorted order of normalized_heatmap_data
filtered_heatmap_data_raw_for_annot = filtered_heatmap_data_raw_for_annot.loc[normalized_heatmap_data.index, :]


print(f"  - Filtered and Normalized heatmap_data DataFrame shape: {normalized_heatmap_data.shape}")
print(f"  - Head of normalized_heatmap_data:")
print(normalized_heatmap_data.head())
print(f"  - Number of NaN values in normalized_heatmap_data: {normalized_heatmap_data.isnull().sum().sum()}")


# --- PLOTTING THE HEATMAP ---
print(f"\n--- Generating heatmap with normalized individual differences ---")
figure_height = len(normalized_heatmap_data) * 0.5 + 3 # Adjust height based on number of traits
figure_width = len(normalized_heatmap_data.columns) * 2.0 + 3 # Adjust width based on number of individuals
if figure_height < 6: figure_height = 6
if figure_width < 10: figure_width = 10

plt.figure(figsize=(figure_width, figure_height))

ax = sns.heatmap(
    normalized_heatmap_data, # Colors are based on Z-scores
    cmap="RdBu_r",           # Red-blue diverging colormap (reversed for high=red, low=blue)
    annot=False,             # <--- SET TO FALSE TO REMOVE NUMBERS
    # fmt=".0f",             # <-- REMOVE THIS LINE
    linewidths=.5,           # Lines between cells
    center=0,                # Center the colormap at 0
    linecolor='black',       # Color of lines between cells
    cbar_kws={'label': 'Normalized Score Difference (Z-score per individual)'}, # Colorbar label
    # annot_kws={"fontsize": 15} # <-- REMOVE THIS LINE
)

# Now, get the colorbar object and set its label's fontweight
cbar = ax.collections[0].colorbar
# Use the direct label string, as 'get_label()' does not exist on the Colorbar object itself
cbar.set_label('Normalized Score Difference (Z-score per individual)', rotation=270, labelpad=15, fontweight='bold')

plt.title(
    f"Heatmap of PGS Score Differences per Individual\n(Colors: Normalized Z-score; Text: Raw Difference - Top {len(normalized_heatmap_data)} Traits by {PERCENTILE_THRESHOLD}th Percentile Absolute Change)",
    fontsize=14, fontweight='bold')
plt.xlabel("Individual ID", fontsize=18, fontweight='bold')
plt.ylabel("PGS Trait", fontsize=18, fontweight='bold')
plt.xticks(fontsize=16, rotation=45, ha='right', fontweight='bold') # Changed ha='center' to 'right' for better rotation alignment
plt.yticks(fontsize=10, fontweight='bold')

plt.tight_layout()

plt.savefig(output_heatmap_filename, format='png', dpi=300, bbox_inches='tight')
print(f"  - Heatmap saved to '{output_heatmap_filename}'")
plt.show()
plt.close()
print(f"--- Plotting complete ---")
