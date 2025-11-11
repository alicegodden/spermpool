# Title: Average PGS across donors z-score normalised
# Author: Dr. Alice M. Godden

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
#input_files = [
#    "PGSscores_difference_D1T4-D1T0.txt",
#    "PGSscores_difference_D2T4-D2T0.txt",
#     "PGSscores_difference_D4T4-D4T0.txt",
#    "PGSscores_difference_D6aT4-D6T0a.txt",
#    "PGSscores_difference_D6bT4-D6bT0.txt",
#]

# MC (active)
input_files = [
    "PGSscores_difference_MC_M8_C-O.txt",
    "PGSscores_difference_MC_M11_C-O.txt",
   "PGSscores_difference_MC_M12_C-O.txt",
    "PGSscores_difference_MC_M13_C-O.txt",
]

phenocodes_path = "phenocodes"
output_heatmap_filename = "PGSscore_heatmap_Average_Donor_Zscore_20_MC_mean.png" # Updated output filename
PERCENTILE_THRESHOLD = 99
MAX_TRAITS_TO_PLOT = 20

# --- LOAD PHENOCODES (Skipping for brevity, assume success) ---
print(f"\n--- Loading phenocodes from: {phenocodes_path} ---")
try:
    phenos = pd.read_csv(phenocodes_path, sep='\t')
    if 'filename' not in phenos.columns or 'description' not in phenos.columns:
        print(f"Error: 'filename' or 'description' column missing in '{phenocodes_path}'.")
        exit(1)
except FileNotFoundError:
    print(f"Error: Phenocodes file not found at '{phenocodes_path}'.")
    exit(1)
except Exception as e:
    print(f"Error loading phenocodes: {e}")
    exit(1)

# --- PROCESS ALL INDIVIDUAL FILES (Skipping for brevity, assume success) ---
all_individuals_dfs = []
print(f"\n--- Processing specified individual files ---")

if not input_files:
    print(f"Error: No input files specified in the 'input_files' list. Exiting.")
    exit(1)

for filepath in input_files:
    if os.path.exists(filepath):
        processed_df = process_individual_file(filepath, phenos)
        if processed_df is not None:
            all_individuals_dfs.append(processed_df)
    else:
        print(f"Warning: Input file not found: '{filepath}'. Skipping.")

if not all_individuals_dfs:
    print("No valid individual data found from input files. Exiting.")
    exit(0)

# Combine all individual DataFrames into a single raw DataFrame
heatmap_data_raw = pd.concat(all_individuals_dfs, axis=1)

# --- TRAIT SELECTION (Filtering based on overall max absolute difference) ---
print(f"\n--- Trait Selection (Percentile Filtering) ---")
max_abs_diff_per_trait = heatmap_data_raw.replace([np.inf, -np.inf], np.nan).stack().abs().groupby(level=0).max().sort_values(ascending=False)
max_abs_diff_per_trait = max_abs_diff_per_trait.dropna()

total_unique_traits = len(max_abs_diff_per_trait)

if total_unique_traits > 0:
    percentile_value = max_abs_diff_per_trait.quantile(1 - (PERCENTILE_THRESHOLD / 100.0))
    if percentile_value == 0.0:
        qualifying_traits_series = max_abs_diff_per_trait[max_abs_diff_per_trait > 0]
    else:
        qualifying_traits_series = max_abs_diff_per_trait[max_abs_diff_per_trait >= percentile_value]

    selected_labels = qualifying_traits_series.head(MAX_TRAITS_TO_PLOT).index.tolist()
else:
    selected_labels = []

if not selected_labels:
    print("No traits met the selection criteria for plotting. Cannot create heatmap. Exiting.")
    exit(0)

# Filter the raw individual data to include only these selected labels
filtered_raw_data = heatmap_data_raw.loc[selected_labels].copy()


# ----------------------------------------------------------------------
# --- NEW CORE LOGIC: CALCULATE, NORMALIZE, AND PREPARE AVERAGE COLUMN ---
# ----------------------------------------------------------------------
print(f"\n--- Calculating and Normalizing Donor Average Z-score ---")

# 1. Calculate the mean across all individual columns for each trait (row-wise)
# This results in the raw average score difference for each phenotype.
donor_average_raw = filtered_raw_data.mean(axis=1).to_frame(name="Donor_Average_Raw")

# 2. Z-SCORE NORMALIZATION (Column-wise, which is across all traits)
def zscore_normalize(series):
    # Skips NaNs automatically
    if series.std() == 0:
        return series * 0.0
    else:
        return (series - series.mean()) / series.std()

# This is the final data for plotting: a single column of Z-scores
normalized_heatmap_data = donor_average_raw.apply(zscore_normalize, axis=0).rename(
    columns={'Donor_Average_Raw': 'Average_Donor_Zscore'}
)

# 3. Final data cleanup and sorting
normalized_heatmap_data.dropna(axis=0, how='all', inplace=True)

if normalized_heatmap_data.empty:
    print("No data remaining after normalization. Cannot create heatmap. Exiting.")
    exit(0)

# Sort traits by their normalized average score
normalized_heatmap_data = normalized_heatmap_data.sort_values(
    by='Average_Donor_Zscore', ascending=False)


print(f"  - Final heatmap data shape: {normalized_heatmap_data.shape}")
print(f"  - Head of final normalized data:")
print(normalized_heatmap_data.head())


# --- PLOTTING THE HEATMAP (Single Column) ---
print(f"\n--- Generating single-column heatmap (Average Z-score) ---")

# Use a fixed, narrow figure width since there is only one column.
# Height is dynamically set based on the number of traits for good readability.
figure_width = 4.0
figure_height = len(normalized_heatmap_data) * 0.5 + 2.5 # Adjusted for traits (rows)
if figure_height < 6: figure_height = 6

plt.figure(figsize=(figure_width, figure_height))

ax = sns.heatmap(
    normalized_heatmap_data,
    cmap="RdBu_r",
    annot=False,
    linewidths=.5,
    center=0,
    linecolor='black',
    cbar_kws={'label': 'Normalized Average Score Difference (Z-score)'},
)

# Set Colorbar label fontweight
cbar = ax.collections[0].colorbar
cbar.set_label('Normalized Average Score Difference (Z-score)', rotation=270, labelpad=15, fontweight='bold')

# Only one column, so set the x-tick label manually
ax.set_xticks([0.5])
ax.set_xticklabels(['Donor Average'], fontsize=12, rotation=45, ha='right', fontweight='bold')
ax.setyticklabels(fontsize=15, rotation=30, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')


# Update Title
plt.title(
    f"Average PGS Score Z-score (Across Donors)\n(Top {len(normalized_heatmap_data)} Traits)",
    fontsize=20, fontweight='bold')
plt.xlabel("", fontsize=18, fontweight='bold') # Remove x-label, as the tick label is sufficient
plt.ylabel("PGS Trait", fontsize=18, fontweight='bold')

plt.tight_layout()

plt.savefig(output_heatmap_filename, format='png', dpi=300, bbox_inches='tight')
print(f"  - Heatmap saved to '{output_heatmap_filename}'")
plt.show()
plt.close()
print(f"--- Plotting complete ---")
