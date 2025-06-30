# matching rows from your vcf files to PGS file, multiplying the AAF with the PGS scores to generate one final score
import os
import gzip
import csv
import glob
from multiprocessing import Pool, cpu_count
from tqdm import tqdm  #   Progress bar

# === CONFIG ===
aaf_file = "/allele_freq/GRCh37_sample.aaf.tsv"
input_dir = "filtered"
output_file = "GRCh37_sample_all_scores.tsv.gz"

# === Load AAF file into memory ===
def load_aaf_dict(path):
    aaf_dict = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row['CHROM'], row['POS'], row['REF'], row['ALT'])
            try:
                aaf_dict[key] = float(row['AAF'])
            except ValueError:
                continue
    return aaf_dict

aaf_variants = load_aaf_dict(aaf_file)

# === Function to score one file ===
def process_file(filepath):
    product_sum = 0.0
    match_file = os.path.basename(filepath)

    try:
        with gzip.open(filepath, 'rt') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for row in reader:
                key = (row.get('chr'), row.get('pos'), row.get('ref'), row.get('alt'))
                try:
                    beta = float(row.get('beta_EUR'))
                    aaf = aaf_variants.get(key)
                    if aaf is not None:
                        product_sum += aaf * beta
                except (TypeError, ValueError):
                    continue
    except Exception as e:
        print(f"Error in {match_file}: {e}")
        return (match_file, "ERROR")

    return (match_file, product_sum)

# === Main execution ===
if __name__ == "__main__":
    input_files = glob.glob(os.path.join(input_dir, "*.tsv.bgz"))
    num_workers = min(cpu_count(), 24)

    print(f"Processing {len(input_files)} files using {num_workers} workers...")

    with Pool(num_workers) as pool:
        results = list(tqdm(pool.imap_unordered(process_file, input_files), total=len(input_files)))

    print("Writing results...")

    with gzip.open(output_file, 'wt', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['match_file', 'product_total'])
        for match_file, total in results:
            writer.writerow([match_file, total])
