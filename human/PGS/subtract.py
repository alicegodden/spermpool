# Title: Subtracting treatment from control polygenic risk scores
# Author: Dr. Alice M. Godden

import pandas as pd
import gzip

# === Input files ===
file1 = "GRCh37_D6T4a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz" #t4 or outer
file2 = "GRCh37_D6T0a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz" # t0 or centre
output_file = "PGSscores_difference_D6aT4-D6T0a.txt"  # Plain text

# === Load both files ===
def load_gz_tsv(file_path):
    with gzip.open(file_path, 'rt') as f:
        return pd.read_csv(f, sep="\t")

df1 = load_gz_tsv(file1)
df2 = load_gz_tsv(file2)

# === Align by match_file ===
df1 = df1.set_index('match_file')
df2 = df2.set_index('match_file')

# Inner join (only rows that exist in both)
merged = df1.join(df2, lsuffix='_T4', rsuffix='_T0', how='inner')

# === Compute difference ===
merged['score_difference'] = merged['product_total_T4'] - merged['product_total_T0']

# === Output result to plain text file ===
output_df = merged[['score_difference']].reset_index()
output_df.to_csv(output_file, sep='\t', index=False)

print(f"✅ Done! Output saved to: {output_file}")
