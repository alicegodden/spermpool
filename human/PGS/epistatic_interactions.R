import pandas as pd
import numpy as np
import os
import re
from scipy.stats import pearsonr
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

# --- Configuration ---
INPUT_DIR = '.'  # Current directory
PHENOCODES_FILE = 'phenocodes'  # File containing phenotype code mappings
# List of input score files
INPUT_FILES = [
    'GRCh37_D1T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D1T4_merged.dedup_RG_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D2T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D2T4_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D4T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D4T4_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D6T0a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D6T4a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D6T0b_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz',
    'GRCh37_D6T4b_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz'
]

# Regex to extract SampleID from filename, e.g., 'D1T0', 'D1T4', 'D6T0a', 'D6T4b'
# This pattern is more robust to handle D#T# and D#T#a/b formats
sample_id_pattern = re.compile(r'GRCh37_(D\d+T[04][ab]?)_merged')

OUTPUT_CORR_MATRIX_R = 'pearson_r_matrix_change_scores.tsv'
OUTPUT_CORR_MATRIX_P = 'pearson_p_matrix_change_scores.tsv'
OUTPUT_SIGNIFICANT_PAIRS = 'significant_correlations_change_scores_SU.tsv'
OUTPUT_HEATMAP_R = 'correlation_heatmap_R_su.png'
OUTPUT_HEATMAP_P = 'correlation_heatmap_P_su.png'

P_VALUE_THRESHOLD = 0.05
FDR_THRESHOLD = 0.05  # For False Discovery Rate correction
TOP_PHENOTYPES_COUNT = 30  # Number of top phenotypes to select for global correlation


# --- Functions ---

def load_phenocodes(filepath):
    """Loads phenocode mappings from a file."""
    try:
        phenocodes_df = pd.read_csv(filepath, sep='\t', header=None, names=['Phenotype_Code', 'Phenotype_Label'])
        print(f"  - Successfully loaded phenocodes. Rows: {phenocodes_df.shape[0]}, Cols: {phenocodes_df.shape[1]}")
        return phenocodes_df
    except FileNotFoundError:
        print(f"Error: Phenocodes file not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error loading phenocodes file: {e}")
        return None


def extract_sample_id(filename):
    """Extracts SampleID (e.g., D1T0, D2T4, D6T0a) from the filename."""
    match = sample_id_pattern.search(filename)
    if match:
        extracted_id = match.group(1)
        return extracted_id
    else:
        return None


def load_and_process_scores(filepath, phenocodes_df):
    """Loads a single score file, maps phenocodes, and handles duplicates."""
    sample_id = extract_sample_id(os.path.basename(filepath))
    if sample_id is None:
        print(f"  - Warning: Could not extract SampleID from filename '{os.path.basename(filepath)}'. Skipping.")
        return pd.DataFrame()  # Return empty DataFrame if SampleID cannot be extracted

    try:
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
        # Renaming 'match_file' to 'Phenotype_Code' and 'product_total' to 'Score'
        df.rename(columns={'match_file': 'Phenotype_Code', 'product_total': 'Score'}, inplace=True)

        # Merge with phenocodes to get labels
        df = pd.merge(df, phenocodes_df, on='Phenotype_Code', how='left')

        # Handle cases where Phenotype_Code might not have a label
        # FIX: Address FutureWarning by not using inplace=True with chained assignment
        df['Phenotype_Label'] = df['Phenotype_Label'].fillna(df['Phenotype_Code'])

        # Check for and handle duplicate trait labels for the same SampleID
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
    Restructures the data and calculates change scores (T4 - T0) for each donor and phenotype.
    Handles 'D#T#a/b' as distinct donor identifiers for pairing.
    """
    # Create a unique donor identifier that handles 'D1', 'D6a', 'D6b' etc.
    # This extracts 'D1', 'D6a', 'D6b' from 'D1T0', 'D6T0a', 'D6T0b' etc.
    all_raw_data_long['donor_id_for_pairing'] = all_raw_data_long['SampleID'].str.replace('T0', '').str.replace('T4',
                                                                                                                '')
    all_raw_data_long['Timepoint'] = all_raw_data_long['SampleID'].str.extract(r'(T[04][ab]?)')

    # Pivot the table to get T0 and T4 (and a/b variants) scores side-by-side
    df_pivot_revised = all_raw_data_long.pivot_table(
        index=['donor_id_for_pairing', 'Phenotype_Label'],
        columns='Timepoint',
        values='Score'
    ).reset_index()

    merged_scores = pd.DataFrame()

    # List of timepoint pairs to look for
    timepoint_pairs = [('T0', 'T4'), ('T0a', 'T4a'), ('T0b', 'T4b')]

    for t0_col, t4_col in timepoint_pairs:
        if t0_col in df_pivot_revised.columns and t4_col in df_pivot_revised.columns:
            temp_df = df_pivot_revised[['donor_id_for_pairing', 'Phenotype_Label', t0_col, t4_col]].copy()
            temp_df['Change_Score'] = temp_df[t4_col] - temp_df[t0_col]
            temp_df.rename(columns={t0_col: 'T0', t4_col: 'T4'}, inplace=True)  # Standardize column names
            merged_scores = pd.concat([merged_scores, temp_df], ignore_index=True)

    # Ensure df_paired_scores is always defined, even if empty
    df_paired_scores = merged_scores.dropna(subset=['Change_Score'])

    # Ensure 'donor_id' column is present as 'donor_id_for_pairing' for consistency downstream
    df_paired_scores.rename(columns={'donor_id_for_pairing': 'donor_id'}, inplace=True)

    print(f"  - Data restructured. Total paired change score entries: {df_paired_scores.shape[0]}")

    return df_paired_scores


def generate_heatmaps(R_matrix, P_matrix, output_r_path, output_p_path, top_phenotypes):
    """Generates and saves heatmaps for R and P matrices."""

    # Filter matrices to include only the top phenotypes
    R_matrix_filtered = R_matrix.loc[top_phenotypes, top_phenotypes]
    P_matrix_filtered = P_matrix.loc[top_phenotypes, top_phenotypes]

    plt.figure(figsize=(16, 14))
    sns.heatmap(R_matrix_filtered, annot=False, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title('Pearson Correlation (R) Matrix of Change Scores (Top Phenotypes)')
    plt.tight_layout()
    plt.savefig(output_r_path)
    plt.close()
    print(f"  - R-value heatmap saved to {output_r_path}")

    plt.figure(figsize=(16, 14))
    sns.heatmap(P_matrix_filtered, annot=False, cmap='viridis_r', fmt=".2e",
                linewidths=.5)  # Use _r for lower p to be darker
    plt.title('P-value Matrix of Change Scores (Top Phenotypes)')
    plt.tight_layout()
    plt.savefig(output_p_path)
    plt.close()
    print(f"  - P-value heatmap saved to {output_p_path}")


# --- Main Script ---
if __name__ == "__main__":
    print(f"DEBUG: Current working directory: {os.getcwd()}")
    print(f"DEBUG: Input files to load: {', '.join(INPUT_FILES)}")

    # 1. Load Phenocodes
    print("\n--- Loading phenocodes from: phenocodes ---")
    phenocodes_df = load_phenocodes(os.path.join(INPUT_DIR, PHENOCODES_FILE))
    if phenocodes_df is None:
        exit("Exiting due to phenocodes loading error.")

    # 2. Read all raw score data
    print("\n--- Reading all raw score data from input files ---")
    all_raw_data_list = []
    for f in INPUT_FILES:
        full_path = os.path.join(INPUT_DIR, f)
        processed_df = load_and_process_scores(full_path, phenocodes_df)
        if not processed_df.empty:
            all_raw_data_list.append(processed_df)

    if not all_raw_data_list:
        exit("No valid score data loaded. Exiting.")

    all_raw_data_long = pd.concat(all_raw_data_list, ignore_index=True)
    print(f"  - Total raw data rows loaded: {all_raw_data_long.shape[0]}")

    # 3. Calculate Change Scores
    print("\n--- Starting Cross-Phenotype Correlation Analysis of Change Scores ---")
    print("  - Restructuring data to calculate change scores per donor and phenotype...")
    df_paired_scores = calculate_change_scores(all_raw_data_long)

    if df_paired_scores.empty:
        exit("No paired change scores found. Cannot perform correlation. Exiting.")

    # Pivot to wide format for correlation, with donors as rows and phenotypes as columns
    # The 'donor_id' column is now the combined 'D#', 'D#a', 'D#b' from `calculate_change_scores`
    df_change_scores_wide = df_paired_scores.pivot(
        index='donor_id',
        columns='Phenotype_Label',
        values='Change_Score'
    )

    # 4. Select top phenotypes for analysis (e.g., those with most data points or largest variance)
    # For simplicity, selecting phenotypes with the most non-NA change scores across donors
    print("  - Selecting top phenotypes for global correlation analysis...")
    phenotype_completeness = df_change_scores_wide.count().sort_values(ascending=False)

    # Filter to phenotypes that have at least 2 non-NaN values for correlation calculation
    eligible_phenotypes = phenotype_completeness[phenotype_completeness >= 2].index

    if len(eligible_phenotypes) == 0:
        exit("No phenotypes with sufficient data points (at least 2 donors) for correlation. Exiting.")

    # Prioritize eligible phenotypes
    top_phenotypes_to_consider = eligible_phenotypes[:TOP_PHENOTYPES_COUNT].tolist()

    # Filter df_change_scores_wide to only include these top phenotypes
    df_change_scores_wide = df_change_scores_wide[top_phenotypes_to_consider]

    print(f"  - Selected {df_change_scores_wide.shape[1]} top phenotypes for correlation.")

    # 5. Perform Correlation Analysis
    print("  - Calculating Pearson correlation matrix and p-values...")

    # Calculate Pearson R matrix using pandas built-in .corr()
    R_matrix = df_change_scores_wide.corr(method='pearson')

    # Calculate P-value matrix manually using scipy.stats.pearsonr
    pheno_cols = df_change_scores_wide.columns
    P_matrix = pd.DataFrame(np.nan, index=pheno_cols, columns=pheno_cols)  # Initialize P-value matrix with NaN

    # Store all p-values for FDR correction
    all_p_values_for_fdr = []

    # Generate unique pairs for p-value calculation
    unique_pairs = list(itertools.combinations(pheno_cols, 2)) + [(col, col) for col in pheno_cols]

    for col1_name, col2_name in unique_pairs:
        # Extract data for the current pair, ensuring it's a Series
        series1 = df_change_scores_wide[col1_name].dropna()
        series2 = df_change_scores_wide[col2_name].dropna()

        # Find common indices (donors) that have valid data for both series
        common_donors = series1.index.intersection(series2.index)

        if len(common_donors) >= 2:  # Pearsonr requires at least 2 data points with common observations
            data1 = series1.loc[common_donors]
            data2 = series2.loc[common_donors]

            _, p_value = pearsonr(data1, data2)

            # Ensure p_value is a scalar float
            if isinstance(p_value, (np.ndarray, pd.Series)) and p_value.ndim == 0:
                p_value = p_value.item()  # Extract scalar from 0-dim array/series
            elif not isinstance(p_value, (float, np.floating)):
                print(
                    f"ERROR: Unexpected p_value type from pearsonr: {type(p_value)} for {col1_name} vs {col2_name}. Setting to NaN.")
                p_value = np.nan  # Fallback if p_value is not a float or 0-dim numpy array

            P_matrix.loc[col1_name, col2_name] = float(p_value)  # Ensure final assignment is float
            P_matrix.loc[col2_name, col1_name] = float(p_value)  # Fill the symmetrical position

            if col1_name != col2_name:  # Only add unique off-diagonal p-values for FDR
                all_p_values_for_fdr.append(p_value)

        elif col1_name == col2_name:  # For self-correlation (diagonal), p-value is 0.0
            P_matrix.loc[col1_name, col2_name] = 0.0

            # 6. Apply FDR Correction
    print(f"  - Applying FDR correction (Benjamini-Hochberg) with alpha = {FDR_THRESHOLD}...")
    if all_p_values_for_fdr:
        # Note: For simplicity and based on previous discussions, the significant correlations
        # table and heatmaps still use raw p-values. If you need a fully FDR-corrected
        # P_matrix or significance table, the mapping back from pvals_corrected to
        # the specific pairs would need to be implemented.
        reject, pvals_corrected, _, _ = multipletests(all_p_values_for_fdr, alpha=FDR_THRESHOLD, method='fdr_bh')
        # You could use `reject` here to filter significant pairs if needed for a separate output.

    # 7. Report Significant Correlations (based on P_VALUE_THRESHOLD)
    print(f"  - Reporting significant correlations (P < {P_VALUE_THRESHOLD})...")
    significant_correlations = []
    for i in range(len(pheno_cols)):
        for j in range(i + 1, len(pheno_cols)):  # Only check upper triangle to avoid duplicates and self-correlation
            col1_name = pheno_cols[i]
            col2_name = pheno_cols[j]

            p_value = P_matrix.loc[col1_name, col2_name]
            r_value = R_matrix.loc[col1_name, col2_name]

            if not pd.isna(p_value) and p_value < P_VALUE_THRESHOLD:
                significant_correlations.append({
                    'Phenotype_1': col1_name,
                    'Phenotype_2': col2_name,
                    'Pearson_R': r_value,
                    'P_Value': p_value
                })

    df_significant = pd.DataFrame(significant_correlations)
    if not df_significant.empty:
        df_significant = df_significant.sort_values(by='P_Value')
        df_significant.to_csv(os.path.join(INPUT_DIR, OUTPUT_SIGNIFICANT_PAIRS), sep='\t', index=False)
        print(f"  - Significant correlations saved to {os.path.join(INPUT_DIR, OUTPUT_SIGNIFICANT_PAIRS)}")
    else:
        print("  - No significant correlations found at the specified P-value threshold.")

    # 8. Save Correlation Matrices
    R_matrix.to_csv(os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_R), sep='\t')
    P_matrix.to_csv(os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_P), sep='\t')
    print(f"  - Pearson R matrix saved to {os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_R)}")
    print(f"  - P-value matrix saved to {os.path.join(INPUT_DIR, OUTPUT_CORR_MATRIX_P)}")

    # 9. Generate Heatmaps
    print("\n--- Generating Heatmaps ---")
    if len(top_phenotypes_to_consider) > 1:  # Need at least 2 phenotypes for a meaningful heatmap
        generate_heatmaps(R_matrix, P_matrix,
                          os.path.join(INPUT_DIR, OUTPUT_HEATMAP_R),
                          os.path.join(INPUT_DIR, OUTPUT_HEATMAP_P),
                          top_phenotypes_to_consider)
    else:
        print("  - Not enough phenotypes selected to generate heatmaps (need at least 2).")

    print("\n--- Script finished ---")



                                                                  ##PLOTTING
                                                                  # manually relabel phenotypes to meaningful label
                                                                  import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os  # Import os module for path manipulation

# --- Configuration (copied from your script) ---
TOP_PHENOTYPES_COUNT = 30  # Number of top phenotypes to select for heatmap


# --- Functions ---

def load_correlation_data(filepath, separator='\t'):
    """
    Loads the correlation data from a TSV/CSV file.
    Assumes the file contains 'Phenotype_1_Description', 'Phenotype_2_Description',
    'Pearson_R', and 'P_Value' columns.
    """
    try:
        df = pd.read_csv(filepath, sep=separator)
        required_cols = ['Phenotype_1_Description', 'Phenotype_2_Description', 'Pearson_R', 'P_Value']

        # Check if all required columns are present
        if not all(col in df.columns for col in required_cols):
            missing_cols = [col for col in required_cols if col not in df.columns]
            raise ValueError(f"Input file is missing required columns: {missing_cols}\n"
                             f"Please ensure your file has all of: {required_cols}.")

        # Ensure numeric types for Pearson_R and P_Value, coercing errors to NaN
        df['Pearson_R'] = pd.to_numeric(df['Pearson_R'], errors='coerce')
        df['P_Value'] = pd.to_numeric(df['P_Value'], errors='coerce')

        # Drop rows where critical numeric conversions failed (i.e., 'Pearson_R' or 'P_Value' are NaN)
        df.dropna(subset=['Pearson_R', 'P_Value'], inplace=True)

        return df
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{filepath}' is empty.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading the file: {e}")
        return None


def prepare_matrix(df, value_column):
    """
    Transforms the long-format correlation or p-value data (pairs) into a square matrix.
    """
    # Get all unique phenotype descriptions to form the matrix axes
    all_phenotypes = pd.concat([df['Phenotype_1_Description'], df['Phenotype_2_Description']]).unique()

    # Initialize an empty square DataFrame with NaN values
    matrix = pd.DataFrame(index=all_phenotypes, columns=all_phenotypes, dtype=float)

    # Populate the matrix with values, ensuring symmetry
    for _, row in df.iterrows():
        pheno1_desc = row['Phenotype_1_Description']
        pheno2_desc = row['Phenotype_2_Description']
        value = row[value_column]

        matrix.loc[pheno1_desc, pheno2_desc] = value
        matrix.loc[pheno2_desc, pheno1_desc] = value  # Populate the symmetric entry

    return matrix


def plot_heatmap(matrix, title, output_path, annot, cmap, fmt, linewidths):
    """
    Generates and saves a heatmap of the provided matrix, matching the user's style.
    """
    plt.figure(figsize=(16, 14))  # Use user-specified figure size

    # Create a mask for the upper triangle if the matrix is square
    # This avoids redundant information in a symmetric correlation matrix
    mask = None
    if matrix.shape[0] == matrix.shape[1]:  # Check if the matrix is square
        mask = np.triu(np.ones_like(matrix, dtype=bool))

    # Generate the heatmap using seaborn
    sns.heatmap(
        matrix,
        mask=mask,  # Apply the mask
        annot=annot,  # Use provided annotation setting (False as per your snippet)
        fmt=fmt,  # Use provided format string (e.g., ".2f", ".2e")
        cmap=cmap,  # Use the specified colormap (e.g., 'coolwarm', 'viridis_r')
        linewidths=linewidths  # Use the specified linewidths (0.5 as per your snippet)
        # vmin/vmax are not specified in your snippet, so we omit them for exact style matching
    )

    plt.title(title, fontsize=16)  # Set title as per your snippet
    plt.xticks(rotation=90, ha='right')  # Rotate x-axis labels for readability
    plt.yticks(rotation=0)  # Keep y-axis labels horizontal
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping, as per your snippet

    # Save the plot to the specified output path and close it
    plt.savefig(output_path)
    plt.close()
    print(f"  - Heatmap saved to {output_path}")


# --- Main Script ---
if __name__ == "__main__":
    print("--- Pearson Correlation and P-value Heatmap Generator (Style-Matched) ---")
    print("This script requires your input file to have the following columns:")
    print("  'Phenotype_1_Description'")
    print("  'Phenotype_2_Description'")
    print("  'Pearson_R'")
    print("  'P_Value'")

    # Prompt user for the input file path
    input_filepath = input(
        "\nEnter the path to your correlation data file (e.g., 'my_correlations_with_descriptions.tsv'): ")

    # Prompt user for the file separator
    file_separator = input("Enter the file separator (e.g., '\\t' for TSV, ',' for CSV): ")
    if file_separator == '\\t':  # Handle the literal string '\t'
        file_separator = '\t'

    # Define output directory and file paths for the heatmaps
    output_dir = "heatmaps_styled"  # Creates a new directory for these styled heatmaps
    os.makedirs(output_dir, exist_ok=True)  # Create the output directory if it doesn't exist
    output_r_path = os.path.join(output_dir, "correlation_heatmap_R.png")  # Matches your output filename
    output_p_path = os.path.join(output_dir, "correlation_heatmap_P.png")  # Matches your output filename

    # Load the data
    df = load_correlation_data(input_filepath, separator=file_separator)

    if df is not None:
        if df.empty:
            print("The loaded file is empty or contains no valid data rows after processing. Cannot generate heatmaps.")
        else:
            all_phenotypes_in_data = pd.concat([df['Phenotype_1_Description'], df['Phenotype_2_Description']]).unique()
            if len(all_phenotypes_in_data) < 2:
                print(
                    f"Warning: Only {len(all_phenotypes_in_data)} unique phenotype(s) found in your data. A heatmap requires at least two distinct phenotypes to show relationships.")
                print("Please check your input data for sufficient unique phenotype pairs.")
            else:
                # Prepare Pearson R-value matrix
                print("\nPreparing Pearson R-value matrix...")
                R_matrix = prepare_matrix(df, 'Pearson_R')
                # For R-value matrix, self-correlation is typically 1.0. Fill diagonal if not present in input.
                np.fill_diagonal(R_matrix.values, 1.0)

                # Prepare P-value matrix
                print("Preparing P-value matrix...")
                P_matrix = prepare_matrix(df, 'P_Value')
                # For P-value matrix, diagonal is usually NaN or 0, not 1.0. It's left as is if not in input data.

                # --- Simulate filtering for top phenotypes as in your provided script ---
                # This step selects phenotypes that have at least 2 non-NaN values
                # and then takes the top N by completeness, mirroring your script's logic.
                print("Selecting top phenotypes for heatmap generation...")

                phenotype_completeness = R_matrix.count().sort_values(ascending=False)
                eligible_phenotypes = phenotype_completeness[phenotype_completeness >= 2].index

                # Select the top N phenotypes from the eligible ones
                top_phenotypes_to_consider = eligible_phenotypes[:TOP_PHENOTYPES_COUNT].tolist()

                if len(top_phenotypes_to_consider) < 2:
                    print(
                        f"Only {len(top_phenotypes_to_consider)} top phenotypes selected after filtering. Need at least 2 for heatmap generation. Skipping heatmaps.")
                else:
                    # Filter R_matrix and P_matrix to include only these top phenotypes
                    R_matrix_filtered = R_matrix.loc[top_phenotypes_to_consider, top_phenotypes_to_consider]
                    P_matrix_filtered = P_matrix.loc[top_phenotypes_to_consider, top_phenotypes_to_consider]

                    print(f"  - Selected {len(top_phenotypes_to_consider)} top phenotypes for heatmap generation.")

                    # Plot and save R-value heatmap, matching your exact style
                    print("Generating Pearson R-value heatmap...")
                    plot_heatmap(R_matrix_filtered,
                                 'Pearson Correlation (R) Matrix of Change Scores (Top Phenotypes)',
                                 output_r_path,
                                 annot=False,  # As per your snippet
                                 cmap='coolwarm',  # As per your snippet
                                 fmt=".2f",  # As per your snippet
                                 linewidths=.5)  # As per your snippet

                    # Plot and save P-value heatmap, matching your exact style
                    print("Generating P-value heatmap...")
                    plot_heatmap(P_matrix_filtered,
                                 'P-value Matrix of Change Scores (Top Phenotypes)',
                                 output_p_path,
                                 annot=False,  # As per your snippet
                                 cmap='viridis_r',  # As per your snippet
                                 fmt=".2e",  # As per your snippet
                                 linewidths=.5)  # As per your snippet

                    print("\nHeatmap generation complete. Check the 'heatmaps_styled' directory for your images.")
