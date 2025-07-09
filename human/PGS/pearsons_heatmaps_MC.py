# Title: Pearsons' correlation with FDR correction MethylCellulose
# Author: Dr. Alice M. Godden

import pandas as pd
import numpy as np
import os
import re
from scipy.stats import pearsonr
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import matplotlib.colors

# --- SET INTERACTIVE PLOTTING MODE (optional) ---
plt.ion()

# --- Configuration ---
INPUT_DIR = '.'  # Current directory
PHENOCODES_FILE = 'phenocodes.tsv'  # File containing phenotype code mappings
# List of input score files for M#_C#/O# convention
INPUT_FILES = [
    "GRCh37_M8_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M8_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M11_C2._merged_dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",  # Note: C2 for M11
    "GRCh37_M11_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M12_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M12_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M13_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
    "GRCh37_M13_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
]

# Regex to extract SampleID from filename based on M#_C#/O# pattern
# This refined regex accounts for both '.dedup' and '._merged_dedup' variations
sample_id_pattern = re.compile(r'GRCh37_([MO]\d+_[CO]\d+)(?:\.?_merged_dedup|\.dedup)?_bcftools')

OUTPUT_CORR_MATRIX_R = 'pearson_r_matrix_change_scores_MC.tsv'
OUTPUT_CORR_MATRIX_P = 'pearson_p_matrix_change_scores_MC.tsv'
OUTPUT_SIGNIFICANT_PAIRS = 'significant_correlations_change_scores_MC.tsv'
OUTPUT_HEATMAP_R = 'correlation_heatmap_R_top100_MC.png'
OUTPUT_HEATMAP_P = 'correlation_heatmap_P_top100_MC.png'
OUTPUT_SIGNIFICANT_HEATMAP_R = 'correlation_heatmap_significant_R_MC.png'
OUTPUT_TOP_CORRELATIONS_PLOT = 'top_20_significant_correlations_MC.png'

P_VALUE_THRESHOLD = 0.05
FDR_THRESHOLD = 0.05
TOP_PHENOTYPES_COUNT = 100  # Changed from 75 to 100 for broader initial analysis
NUM_TOP_CORRELATIONS_TO_PLOT = 20


# --- Functions ---

def extract_phenotype_key(filename):
    """
    Extracts a consistent key from the phenotype filename by intelligently
    removing '_filtered.tsv.bgz' or '.tsv.bgz' extensions.
    """
    if filename.endswith('_filtered.tsv.bgz'):
        key = filename.replace('_filtered.tsv.bgz', '')
    elif filename.endswith('.tsv.bgz'):
        key = filename.replace('.tsv.bgz', '')
    else:
        key = filename
    return key


def load_phenocodes(filepath):
    """
    Loads phenocode mappings from a file into a dictionary.
    """
    phenocode_map = {}
    print(f"  - Loading phenocodes from {filepath}")
    try:
        df_pheno = pd.read_csv(filepath, sep='\t', header=0)

        if 'description' not in df_pheno.columns or 'filename' not in df_pheno.columns:
            print(f"Error: phenocodes file '{filepath}' must contain 'description' and 'filename' columns.")
            return None

        for index, row in df_pheno.iterrows():
            description = str(row['description']).strip()
            full_filename = str(row['filename']).strip()

            key = extract_phenotype_key(full_filename)

            if not key:
                print(
                    f"Warning: Row {index + 2}: Could not extract valid key from phenocode filename: '{full_filename}'. Skipping this entry.")
                continue

            if not description:
                print(
                    f"Warning: Row {index + 2}: Empty description for filename '{full_filename}'. Skipping this entry.")
                continue

            phenocode_map[key] = description

            if index < 5:  # Reduced DEBUG output for phenocodes
                print(f"DEBUG: Phenocode Map Sample (Line {index + 2}): Key='{key}', Description='{description}'")

        print(f"  - Successfully loaded {len(phenocode_map)} phenocode descriptions.")
        return phenocode_map
    except FileNotFoundError:
        print(f"Error: Phenocodes file not found at {filepath}")
        return None
    except Exception as e:
        print(f"An error occurred while loading the phenocodes file with pandas: {e}")
        return None


def extract_sample_id(filename):
    """
    Extracts SampleID (e.g., M8_C1, M11_O1) from the filename using the updated regex.
    """
    match = sample_id_pattern.search(filename)
    if match:
        extracted_id = match.group(1)  # This group captures the full M#_C# or M#_O#
        return extracted_id
    else:
        print(f"DEBUG: No SampleID match for file: {filename}")
        return None


def load_and_process_scores(filepath, phenocode_descriptions_map):
    """
    Loads a single score file, maps phenotype filenames to their descriptions,
    and handles duplicates.
    """
    sample_id = extract_sample_id(os.path.basename(filepath))
    if sample_id is None:
        print(f"  - Warning: Could not extract SampleID from filename '{os.path.basename(filepath)}'. Skipping.")
        return pd.DataFrame()

    try:
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
        df.rename(columns={'match_file': 'Phenotype_Code', 'product_total': 'Score'}, inplace=True)

        df['Phenotype_Key_For_Lookup'] = df['Phenotype_Code'].apply(extract_phenotype_key)

        # DEBUG: Print only a few head/tail for brevity
        if not df.empty:
            print(f"DEBUG: Sample {sample_id} - First 2 Phenotype Keys from score file:")
            print(df['Phenotype_Key_For_Lookup'].head(2).tolist())
            # Removed redundant tail output from previous versions for brevity

        df['Phenotype_Label'] = df['Phenotype_Key_For_Lookup'].map(phenocode_descriptions_map)

        unmapped_phenotypes = df[df['Phenotype_Label'].isna()]['Phenotype_Key_For_Lookup'].unique()
        if len(unmapped_phenotypes) > 0:
            print(f"DEBUG: Sample {sample_id} - First 2 UNMAPPED Phenotype Keys (will show as filenames):")
            print(unmapped_phenotypes[:2].tolist())
            print(f"DEBUG: Total unmapped phenotypes in {os.path.basename(filepath)}: {len(unmapped_phenotypes)}")
        # Removed else: print for less verbose output when all mapped

        df['Phenotype_Label'] = df['Phenotype_Label'].fillna(df['Phenotype_Code'])

        df.drop(columns=['Phenotype_Key_For_Lookup'], inplace=True)

        duplicate_labels = df.duplicated(subset=['Phenotype_Label'], keep=False)
        if duplicate_labels.any():
            print(
                f"  - Warning: Duplicate trait labels found in {os.path.basename(filepath)}. Taking mean score for duplicates.")
            df = df.groupby('Phenotype_Label')['Score'].mean().reset_index()

        df['SampleID'] = sample_id
        return df
    except Exception as e:
        print(f"Error processing file {os.path.basename(filepath)}: {e}")
        return pd.DataFrame()


def calculate_change_scores(all_raw_data_long):
    """
    Restructures the data and calculates change scores (O1 - C1/C2) for each donor and phenotype.
    Handles 'M#_C#' and 'M#_O#' as timepoints.
    """
    # DEBUG: Check initial state of all_raw_data_long
    print("\n--- Calculating change scores (O1 - C1/C2) ---")
    print(f"DEBUG: all_raw_data_long shape before processing: {all_raw_data_long.shape}")
    print(f"DEBUG: all_raw_data_long columns: {all_raw_data_long.columns.tolist()}")
    print(f"DEBUG: Unique SampleIDs in all_raw_data_long: {all_raw_data_long['SampleID'].unique().tolist()}")

    # Extract just the donor base ID (e.g., 'M8', 'M11')
    all_raw_data_long['donor_base_id'] = all_raw_data_long['SampleID'].apply(lambda x: x.split('_')[0])

    # Extract the timepoint part (e.g., 'C1', 'C2', 'O1')
    all_raw_data_long['Timepoint'] = all_raw_data_long['SampleID'].apply(lambda x: x.split('_')[1])

    # Pivot the table to get C1/C2 and O1 scores side-by-side
    df_pivot_revised = all_raw_data_long.pivot_table(
        index=['donor_base_id', 'Phenotype_Label'],
        columns='Timepoint',
        values='Score'
    ).reset_index()

    # DEBUG: Check the pivot table structure
    print(f"DEBUG: df_pivot_revised shape: {df_pivot_revised.shape}")
    print(f"DEBUG: df_pivot_revised columns: {df_pivot_revised.columns.tolist()}")
    print(f"DEBUG: df_pivot_revised head:\n{df_pivot_revised.head()}")

    merged_scores = pd.DataFrame()

    # Define the baseline and follow-up timepoints
    baseline_timepoints = ['C1', 'C2']  # C1 or C2 can be baseline
    followup_timepoint = 'O1'

    # Iterate through each donor base ID
    for donor_base in df_pivot_revised['donor_base_id'].unique():
        donor_df = df_pivot_revised[df_pivot_revised['donor_base_id'] == donor_base]

        # Check if an O1 score column exists AND has non-null data for this donor
        if followup_timepoint in donor_df.columns and not donor_df[followup_timepoint].isnull().all():

            # Find the corresponding baseline score (C1 or C2) that exists AND has non-null data
            c_score_col = None  # Store the column name (e.g., 'C1' or 'C2')
            for baseline_col in baseline_timepoints:
                if baseline_col in donor_df.columns and not donor_df[baseline_col].isnull().all():
                    c_score_col = baseline_col
                    break  # Found a baseline, use it

            if c_score_col is not None:
                # Select only relevant columns for this donor, keeping Phenotype_Label as the index
                # This ensures alignment for subtraction
                donor_scores_for_pairing = donor_df[['Phenotype_Label', c_score_col, followup_timepoint]].set_index(
                    'Phenotype_Label')

                # Calculate the change score directly. Pandas aligns by index (Phenotype_Label).
                change_scores = donor_scores_for_pairing[followup_timepoint] - donor_scores_for_pairing[c_score_col]

                # Create the temp_df to concatenate, only including non-NaN change scores
                temp_df_for_concat = pd.DataFrame({
                    'donor_base_id': donor_base,
                    'Phenotype_Label': change_scores.index,
                    'Change_Score': change_scores.values,  # Use .values to avoid index alignment issues during concat
                    'C_Score': donor_scores_for_pairing[c_score_col].values,
                    'O_Score': donor_scores_for_pairing[followup_timepoint].values,
                    'Baseline_Timepoint_Used': c_score_col
                })

                # Only concatenate rows where Change_Score is not NaN
                temp_df_for_concat = temp_df_for_concat.dropna(subset=['Change_Score'])

                if not temp_df_for_concat.empty:
                    merged_scores = pd.concat([merged_scores, temp_df_for_concat], ignore_index=True)
                else:
                    print(f"  - Note: No paired change scores generated for donor {donor_base} after NaN removal.")
            else:
                print(
                    f"  - Warning: No non-null baseline (C1 or C2) score found for donor {donor_base} with O1 scores. Skipping pairing for this donor.")
        else:
            print(f"  - Warning: No non-null O1 score found for donor {donor_base}. Skipping pairing for this donor.")

    # IMPORTANT CHECK BEFORE DROPNAN (already handled by dropna(subset=['Change_Score']) in loop)
    if merged_scores.empty:
        print("  - Resulting 'merged_scores' DataFrame is empty. No change scores calculated.")
        return pd.DataFrame()  # Return empty DataFrame if no data

    # Rename for consistency, handling SettingWithCopyWarning
    df_paired_scores = merged_scores.rename(columns={'donor_base_id': 'donor_id'})

    print(f"  - Data restructured. Total paired change score entries: {df_paired_scores.shape[0]}")
    return df_paired_scores


def generate_heatmaps(R_matrix, P_matrix, output_r_path, output_p_path, top_phenotypes, top_phenotypes_count):
    """Generates and saves heatmaps for R and P matrices using clustermap for clustering."""

    R_matrix_filtered = R_matrix.loc[top_phenotypes, top_phenotypes].copy()
    P_matrix_filtered = P_matrix.loc[top_phenotypes, top_phenotypes].copy()

    # Fill NaN values to ensure clustering can proceed
    # For R-matrix, replace NaNs with 0 (no correlation)
    R_matrix_filtered.fillna(0, inplace=True)
    # For P-matrix, replace NaNs with 1.0 (no significance)
    P_matrix_filtered.fillna(1.0, inplace=True)

    # --- R-value Heatmap (Clustered) ---
    figsize_r_width = max(18, len(top_phenotypes) * 0.2)
    figsize_r_height = max(16, len(top_phenotypes) * 0.2)

    g_r = sns.clustermap(R_matrix_filtered,
                         annot=False,
                         cmap='coolwarm',
                         center=0,
                         vmin=-1, vmax=1,
                         linewidths=.5, linecolor='lightgrey',
                         cbar_kws={'label': 'Pearson Correlation (R)'},
                         figsize=(figsize_r_width, figsize_r_height))

    g_r.fig.suptitle(
        f'Clustered Pearson Correlation (R) Matrix of Change Scores (Top {top_phenotypes_count} Phenotypes)', y=1.02)
    plt.setp(g_r.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(g_r.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    plt.savefig(output_r_path, bbox_inches='tight')
    plt.close()
    print(f"  - Clustered R-value heatmap saved to {output_r_path}")

    # --- P-value Heatmap (Clustered) ---
    figsize_p_width = max(18, len(top_phenotypes) * 0.2)
    figsize_p_height = max(16, len(top_phenotypes) * 0.2)

    g_p = sns.clustermap(P_matrix_filtered,
                         annot=False,
                         cmap='viridis_r',
                         fmt=".2e",
                         linewidths=.5, linecolor='lightgrey',
                         cbar_kws={'label': 'P-value'},
                         figsize=(figsize_p_width, figsize_p_height))

    g_p.fig.suptitle(f'Clustered P-value Matrix of Change Scores (Top {top_phenotypes_count} Phenotypes)', y=1.02)
    plt.setp(g_p.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(g_p.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    plt.savefig(output_p_path, bbox_inches='tight')
    plt.close()
    print(f"  - Clustered P-value heatmap saved to {output_p_path}")


def plot_significant_heatmap(R_matrix, df_significant, output_path, vmin_custom=None, vmax_custom=None,
                             cmap_base_name=None, n_colors_discrete=None):
    """
    Generates a clustered heatmap showing only significantly correlated terms.
    Correlations are filtered by FDR_THRESHOLD.
    This function uses clustermap to automatically group similar phenotypes.
    Accepts vmin_custom and vmax_custom to control the color scale.
    cmap_base_name: Name of a base matplotlib colormap (e.g., 'Reds', 'viridis').
    n_colors_discrete: Number of discrete color steps.
    """
    if df_significant.empty:
        print("  - No significant correlations found to plot heatmap.")
        return

    print(f"  - Generating heatmap for significant correlations...")

    significant_phenotypes = pd.Series(
        df_significant['Phenotype_1'].tolist() + df_significant['Phenotype_2'].tolist()).unique()

    if len(significant_phenotypes) < 2:
        print("  - Not enough unique significant phenotypes (less than 2) to plot heatmap.")
        return

    significant_R_matrix = pd.DataFrame(np.nan, index=significant_phenotypes, columns=significant_phenotypes)

    for _, row in df_significant.iterrows():
        p1 = row['Phenotype_1']
        p2 = row['Phenotype_2']
        r_val = row['Pearson_R']
        significant_R_matrix.loc[p1, p2] = r_val
        significant_R_matrix.loc[p2, p1] = r_val

    np.fill_diagonal(significant_R_matrix.values, 1.0)

    # Fill any remaining NaNs (e.g., from correlations that weren't significant, if any) with 0
    significant_R_matrix.fillna(0, inplace=True)

    annotate_values = len(significant_phenotypes) <= 30

    # Set vmin and vmax for the colormap based on custom inputs or default to -1, 1
    plot_vmin = vmin_custom if vmin_custom is not None else -1
    plot_vmax = vmax_custom if vmax_custom is not None else 1

    # Determine the colormap and its center
    effective_cmap = 'coolwarm'  # Default diverging colormap
    effective_center = 0  # Default center for diverging colormap

    if cmap_base_name and n_colors_discrete:
        # Create a discrete colormap from a base sequential colormap
        base_cmap_continuous = plt.cm.get_cmap(cmap_base_name, 256)  # Get a continuous version
        new_colors = base_cmap_continuous(np.linspace(0, 1, n_colors_discrete))
        effective_cmap = matplotlib.colors.ListedColormap(new_colors)
        effective_center = None  # Center is not meaningful for a sequential discrete map

    g = sns.clustermap(significant_R_matrix,
                       cmap=effective_cmap,  # Use the chosen colormap
                       center=effective_center,  # Use the chosen center
                       vmin=plot_vmin, vmax=plot_vmax,
                       linewidths=.5, linecolor='lightgrey',
                       cbar_kws={'label': 'Pearson Correlation (R)'},
                       annot=annotate_values, fmt=".2f",
                       figsize=(max(15, len(significant_phenotypes) * 0.8), max(12, len(significant_phenotypes) * 0.8)))

    g.fig.suptitle('Significant Pearson Correlations (R) of Change Scores (FDR Corrected)', y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  - Significant correlation heatmap saved to {output_path}")


def plot_top_correlated_pairs(df_significant, num_pairs, output_path):
    """
    Plots the top N most positively and negatively correlated significant phenotype pairs.
    """
    if df_significant.empty:
        print("  - No significant correlations found to plot top pairs.")
        return

    print(f"  - Plotting top {num_pairs} positive and negative significant correlations...")

    df_positive = df_significant[df_significant['Pearson_R'] > 0].sort_values(by='Pearson_R', ascending=False)
    df_negative = df_significant[df_significant['Pearson_R'] < 0].sort_values(by='Pearson_R', ascending=True)

    top_positive = df_positive.head(num_pairs).copy()
    top_negative = df_negative.head(num_pairs).copy()

    if top_positive.empty and top_negative.empty:
        print("  - No positive or negative significant correlations to plot.")
        return

    if not top_positive.empty:
        top_positive['Correlation_Type'] = 'Positive'
    if not top_negative.empty:
        top_negative['Correlation_Type'] = 'Negative'

    plot_df = pd.concat([top_positive, top_negative])

    if plot_df.empty:
        print("  - Combined DataFrame for plotting is empty. No plot generated.")
        return

    plot_df['Pair'] = plot_df['Phenotype_1'] + ' vs. ' + plot_df['Phenotype_2']
    plot_df = plot_df.sort_values(by=['Correlation_Type', 'Pearson_R'], ascending=[False, True])

    plt.figure(figsize=(12, max(6, len(plot_df) * 0.4)))
    sns.barplot(x='Pearson_R', y='Pair', hue='Correlation_Type', data=plot_df,
                palette={'Positive': 'green', 'Negative': 'red'}, dodge=False)

    plt.title(f'Top {num_pairs} Most Significant Positive and Negative Correlations in Change Scores')
    plt.xlabel('Pearson Correlation Coefficient (R)')
    plt.ylabel('Phenotype Pair')
    plt.axvline(x=0, color='grey', linestyle='--', linewidth=0.8)
    plt.legend(title='Correlation Type')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  - Top correlated pairs plot saved to {output_path}")


# --- Main Script ---
if __name__ == "__main__":
    print(f"DEBUG: Current working directory: {os.getcwd()}")
    print(f"DEBUG: Input files to load: {', '.join(INPUT_FILES)}")

    # 1. Load Phenocodes
    print("\n--- Loading phenocodes ---")
    phenocode_descriptions_map = load_phenocodes(os.path.join(INPUT_DIR, PHENOCODES_FILE))
    if phenocode_descriptions_map is None:
        exit("Exiting due to phenocodes loading error. Please check phenocodes.tsv format and presence.")

    # 2. Read all raw score data
    all_raw_data_list = []
    for f in INPUT_FILES:
        full_path = os.path.join(INPUT_DIR, f)
        processed_df = load_and_process_scores(full_path, phenocode_descriptions_map)
        if not processed_df.empty:
            all_raw_data_list.append(processed_df)

    if not all_raw_data_list:
        exit("No valid score data loaded. Exiting.")

    all_raw_data_long = pd.concat(all_raw_data_list, ignore_index=True)

    # 3. Calculate Change Scores
    df_paired_scores = calculate_change_scores(all_raw_data_long)

    if df_paired_scores.empty:
        exit("No paired change scores found. Cannot perform correlation. Exiting.")

    df_change_scores_wide = df_paired_scores.pivot(
        index='donor_id',
        columns='Phenotype_Label',
        values='Change_Score'
    )

    # 4. Select top phenotypes for analysis
    # Filter out phenotypes that don't have at least 2 non-NaN values for correlation
    phenotype_completeness = df_change_scores_wide.count().sort_values(ascending=False)
    eligible_phenotypes = phenotype_completeness[phenotype_completeness >= 2].index

    if len(eligible_phenotypes) == 0:
        exit("No phenotypes with sufficient data points (at least 2 donors) for correlation. Exiting.")

    # Only consider the top N phenotypes for the main R/P matrices and initial heatmaps
    top_phenotypes_to_consider = eligible_phenotypes[:TOP_PHENOTYPES_COUNT].tolist()
    # Ensure all top_phenotypes_to_consider actually exist in the wide dataframe's columns
    # (they should if they came from df_change_scores_wide.columns)
    df_change_scores_wide = df_change_scores_wide[top_phenotypes_to_consider]

    # 5. Perform Correlation Analysis
    print("\n--- Performing correlation analysis ---")
    R_matrix = df_change_scores_wide.corr(method='pearson')

    pheno_cols = df_change_scores_wide.columns
    P_matrix = pd.DataFrame(np.nan, index=pheno_cols, columns=pheno_cols)

    all_raw_p_values_for_fdr_with_pairs = []

    for i in range(len(pheno_cols)):
        for j in range(i, len(pheno_cols)):
            col1_name = pheno_cols[i]
            col2_name = pheno_cols[j]

            series1 = df_change_scores_wide[col1_name].dropna()
            series2 = df_change_scores_wide[col2_name].dropna()
            common_donors = series1.index.intersection(series2.index)

            if col1_name == col2_name:
                P_matrix.loc[col1_name, col2_name] = 0.0
            elif len(common_donors) >= 2:  # Pearson correlation requires at least 2 data points
                data1 = series1.loc[common_donors]
                data2 = series2.loc[common_donors]

                # Handle potential edge case where data becomes constant after subsetting common_donors
                if data1.nunique() < 2 or data2.nunique() < 2:
                    p_value = np.nan  # Cannot calculate correlation if data is constant
                else:
                    _, p_value = pearsonr(data1, data2)

                if isinstance(p_value, (np.ndarray, pd.Series)) and p_value.ndim == 0:
                    p_value = p_value.item()
                elif not isinstance(p_value, (float, np.floating)):
                    p_value = np.nan

                P_matrix.loc[col1_name, col2_name] = float(p_value)
                P_matrix.loc[col2_name, col1_name] = float(p_value)

                # Only add non-NaN p-values for FDR correction
                if not np.isnan(p_value):
                    all_raw_p_values_for_fdr_with_pairs.append({
                        'Phenotype_1': col1_name,
                        'Phenotype_2': col2_name,
                        'Raw_P_Value': p_value
                    })
            else:
                P_matrix.loc[col1_name, col2_name] = np.nan
                P_matrix.loc[col2_name, col1_name] = np.nan

    # 6. Apply FDR Correction
    print("\n--- Applying FDR correction ---")
    df_raw_p_values = pd.DataFrame(all_raw_p_values_for_fdr_with_pairs)
    df_significant_fdr = pd.DataFrame()

    if not df_raw_p_values.empty and 'Raw_P_Value' in df_raw_p_values.columns:
        valid_p_values_series = df_raw_p_values['Raw_P_Value'].dropna()
        valid_indices = valid_p_values_series.index
        valid_p_values = valid_p_values_series.tolist()

        if valid_p_values:
            reject, pvals_corrected, _, _ = multipletests(valid_p_values, alpha=FDR_THRESHOLD, method='fdr_bh')

            df_raw_p_values['FDR_P_Value'] = np.nan
            df_raw_p_values.loc[valid_indices, 'FDR_P_Value'] = pvals_corrected

            df_significant_fdr = df_raw_p_values[df_raw_p_values['FDR_P_Value'] < FDR_THRESHOLD].copy()

            R_lookup = R_matrix.stack().rename('Pearson_R')


            def get_r_value(row):
                try:
                    # Look up R for the specific phenotype pair, handling order
                    # Check (P1, P2) and (P2, P1)
                    if (row['Phenotype_1'], row['Phenotype_2']) in R_lookup.index:
                        return R_lookup.loc[(row['Phenotype_1'], row['Phenotype_2'])]
                    elif (row['Phenotype_2'], row['Phenotype_1']) in R_lookup.index:
                        return R_lookup.loc[(row['Phenotype_2'], row['Phenotype_1'])]
                    else:
                        return np.nan  # Should not happen if phenotypes are from R_matrix
                except KeyError:
                    return np.nan


            df_significant_fdr['Pearson_R'] = df_significant_fdr.apply(get_r_value, axis=1)

            df_significant_fdr.dropna(subset=['Pearson_R'], inplace=True)
            df_significant_fdr.sort_values(by='FDR_P_Value', inplace=True)

            df_significant_fdr.to_csv(os.path.join(INPUT_DIR, OUTPUT_SIGNIFICANT_PAIRS), sep='\t', index=False)
            print(
                f"  - Significant correlations (FDR corrected) saved to {os.path.join(INPUT_DIR, OUTPUT_SIGNIFICANT_PAIRS)}")
        else:
            print("  - No valid p-values to apply FDR correction, or all p-values are NaN.")
    else:
        print("  - No raw p-values found for FDR correction (DataFrame is empty or missing 'Raw_P_Value').")

    # 7. Save Correlation Matrices
    print("\n--- Saving correlation matrices ---")
    R_matrix.to_csv(os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_R), sep='\t')
    P_matrix.to_csv(os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_P), sep='\t')
    print(f"  - Pearson R matrix (all top phenotypes) saved to {os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_R)}")
    print(f"  - P-value matrix (all top phenotypes) saved to {os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_P)}")

    # 8. Generate Heatmaps for all top N phenotypes (CLUSTERED)
    print("\n--- Generating Heatmaps (all top N phenotypes - CLUSTERED) ---")
    if len(top_phenotypes_to_consider) > 1:
        generate_heatmaps(R_matrix, P_matrix,
                          os.path.join(INPUT_DIR, OUTPUT_HEATMAP_R),
                          os.path.join(INPUT_DIR, OUTPUT_HEATMAP_P),
                          top_phenotypes_to_consider,
                          TOP_PHENOTYPES_COUNT)
    else:
        print("  - Not enough phenotypes selected to generate general heatmaps (need at least 2).")

    # 9. Generate Heatmap for only significantly correlated terms
    print("\n--- Generating Heatmap for SIGNIFICANT correlations ---")

    print(
        "WARNING: Custom color scale for significant heatmap (0.95-1.0, 20 steps, Reds) will cause all R < 0.95 to appear as the lowest color.")
    plot_significant_heatmap(R_matrix, df_significant_fdr, os.path.join(INPUT_DIR, OUTPUT_SIGNIFICANT_HEATMAP_R),
                             vmin_custom=-1.0, vmax_custom=1.0,
                             cmap_base_name='coolwarm', n_colors_discrete=1000)

    # 10. Plot Top N significant correlations (using the FDR filtered list)
    print("\n--- Plotting Top Significant Correlations ---")
    plot_top_correlated_pairs(df_significant_fdr, NUM_TOP_CORRELATIONS_TO_PLOT,
                              os.path.join(INPUT_DIR, OUTPUT_TOP_CORRELATIONS_PLOT))

    # --- Print Top N Highest and Lowest Significant R-values to console ---
    if not df_significant_fdr.empty:
        df_positive_top = df_significant_fdr[df_significant_fdr['Pearson_R'] > 0].sort_values(by='Pearson_R',
                                                                                              ascending=False).head(
            NUM_TOP_CORRELATIONS_TO_PLOT)
        df_negative_top = df_significant_fdr[df_significant_fdr['Pearson_R'] < 0].sort_values(by='Pearson_R',
                                                                                              ascending=True).head(
            NUM_TOP_CORRELATIONS_TO_PLOT)

        print(f"\n--- Top {NUM_TOP_CORRELATIONS_TO_PLOT} Highest Significant Pearson R Values ---")
        if not df_positive_top.empty:
            for index, row in df_positive_top.iterrows():
                print(
                    f"  {row['Phenotype_1']} vs. {row['Phenotype_2']}: R = {row['Pearson_R']:.2f}, FDR P = {row['FDR_P_Value']:.2e}")
        else:
            print("  No significant positive correlations found.")

        print(f"\n--- Top {NUM_TOP_CORRELATIONS_TO_PLOT} Lowest Significant Pearson R Values (Most Negative) ---")
        if not df_negative_top.empty:
            for index, row in df_negative_top.iterrows():
                print(
                    f"  {row['Phenotype_1']} vs. {row['Phenotype_2']}: R = {row['Pearson_R']:.2f}, FDR P = {row['FDR_P_Value']:.2e}")
        else:
            print("  No significant negative correlations found.")
    else:
        print("\nNo significant correlations found to list top pairs.")

    print("\n--- Script finished ---")
