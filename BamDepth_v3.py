"""
Purpose: Get the actual depth for varaints in the VCF
Created on Sept 02 23:29:11 2025
Coded by: Vyomesh J
Optimization and testing by: Sophia EA
"""
import os
import gzip
import pandas as pd
import subprocess
import logging
import platform
import getpass
import socket
import psutil
import time
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from io import StringIO
import sys

# Setup logging
logging.basicConfig(
    filename='variant_depth_processing.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Log system information
logging.info("System Information:")
logging.info(f"Hostname: {socket.gethostname()}")
logging.info(f"User: {getpass.getuser()}")
logging.info(f"Platform: {platform.platform()}")
logging.info(f"Processor: {platform.processor()}")
logging.info(f"Python Version: {platform.python_version()}")
logging.info(f"CPU Cores: {psutil.cpu_count(logical=True)}")
logging.info(f"Memory: {round(psutil.virtual_memory().total / (1024**3), 2)} GB")
logging.info(f"Available Memory: {round(psutil.virtual_memory().available / (1024**3), 2)} GB")
logging.info(f"Current CPU Usage: {psutil.cpu_percent(interval=1)}%")

# number of threads
THREADS_PER_SAMPLE = int(os.environ.get("THREADS_PER_SAMPLE", 12))
logging.info(f"--- Starting variant processing with {THREADS_PER_SAMPLE} threads per sample ---")

def process_sample(args):
    sample, Batch_ID, bedfile = args
    sample_temp_dir = os.path.join('tmp', sample)
    os.makedirs(sample_temp_dir, exist_ok=True)

    vcf_path = f"/home/soph/prototype_01/03_Raw_Data_FastQ/OPF001_Test/testdata/cancerPanel/{sample}/{sample}.hard-filtered.vcf"
    bam_path = f"/home/soph/prototype_01/03_Raw_Data_FastQ/OPF001_Test/testdata/cancerPanel/{sample}/{sample}.bam"

    try:
        logging.info(f"Processing sample: {sample}")
        logging.info(f"Current CPU usage: {psutil.cpu_percent(interval=1)}% | Available memory: {round(psutil.virtual_memory().available / (1024**3), 2)} GB")
        start_time = time.time()

        with gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path) as f:
            lines = []
            for line in f:
                if line.startswith('#CHROM'):
                    lines.append(line)
                    break
            for line in f:
                lines.append(line)

# added data type strings
        df_vcf = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=0, dtype=str)
        df_vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample_Info"]
        
# added REF_LEN and ALT_LEN; END calculation
        df_vcf['POS'] = df_vcf['POS'].astype(int)
        df_vcf['REF_LEN'] = df_vcf['REF'].str.len().astype(int)
        df_vcf['ALT_LEN'] = df_vcf['ALT'].str.len().astype(int)
        df_vcf['END'] = df_vcf['POS'] + df_vcf['REF_LEN'] - 1
        df_vcf['End'] = df_vcf['END']

        df_wig_combined = pd.DataFrame()

        def process_variant(row):
            chrom = row['CHROM']
            pos = int(row['POS'])
            ref = row['REF']
            alt = row['ALT']
            ref_len = int(row['REF_LEN'])
            alt_len = int(row['ALT_LEN'])
            end = int(row['END'])

            if ref_len == alt_len:
                # SNP or substitution
                q_start = pos
                q_end = pos
                alt_query_pos = pos
            elif ref_len > alt_len:
                # deletion
                del_start = pos + 1
                del_end = pos + ref_len - 1
                q_start = pos
                q_end = del_end
                alt_query_pos = del_start
            else:
                # insertion
                q_start = pos
                q_end = pos + 1
                alt_query_pos = pos

            pattern = f"{chrom}:{q_start}-{q_end}"
            wig_path = os.path.join(sample_temp_dir, f"{sample}_{q_start}_{q_end}.wig")

            try:
                cmd = [
                    "/home/soph/prototype_01/IGV_2.19.5/igvtools", "count", "-w", "1", "--bases", "--query",
                    pattern, bam_path, wig_path, "hg38"
                ]
                subprocess.run(cmd, check=True)
# check if IGV wig file is present
                if not os.path.exists(wig_path):
                    logging.warning(f"wig file not found for pattern {pattern} (sample {sample})")
                    return None

# reading wig file
                df_w = pd.read_csv(wig_path, sep='\t', skiprows=3, header=None)
                df_w.columns = ['End', 'Base_A', 'Base_C', 'Base_G', 'Base_T', 'Base_N', 'Base_DEL', 'Base_INS']
                df_w['End'] = df_w['End'].astype(int)

                # total depth; excluded Base_DEL to avoid double counting
                df_w['Total_Depth'] = df_w[['Base_A', 'Base_C', 'Base_G', 'Base_T', 'Base_N']].sum(axis=1)

                # find exact coordinate row for alt_query_pos
                target = df_w[df_w['End'] == alt_query_pos]
                if target.empty:
                    # try +/- 1 as a fallback if exact position isn't found
                    target = df_w[df_w['End'].isin([alt_query_pos - 1, alt_query_pos + 1])]
                # extract base counts
                if not target.empty:
                    tr = target.iloc[0]
                    base_a = int(tr['Base_A'])
                    base_c = int(tr['Base_C'])
                    base_g = int(tr['Base_G'])
                    base_t = int(tr['Base_T'])
                    base_n = int(tr['Base_N'])
                    base_del = int(tr['Base_DEL'])
                    base_ins = int(tr['Base_INS'])
                    total_depth = int(tr['Total_Depth'])

                    # alt_count based on variant type (SNP/indel)
                    if ref_len == alt_len:
                        alt_base_col = f"Base_{alt.upper()}"
                        alt_count = int(tr.get(alt_base_col, 0))
                    elif ref_len > alt_len:
                        alt_count = base_del
                    else:
                        alt_count = base_ins
                else:
                    # no exact row found: aggregate across region as fallback
                    base_a = int(df_w['Base_A'].sum())
                    base_c = int(df_w['Base_C'].sum())
                    base_g = int(df_w['Base_G'].sum())
                    base_t = int(df_w['Base_T'].sum())
                    base_n = int(df_w['Base_N'].sum())
                    base_del = int(df_w['Base_DEL'].sum())
                    base_ins = int(df_w['Base_INS'].sum())
                    total_depth = int(df_w['Total_Depth'].sum())

                    if ref_len == alt_len:
                        # SNP fallback: try middle row's base for alt if possible
                        mid_row = df_w.iloc[len(df_w) // 2]
                        alt_count = int(mid_row.get(f"Base_{alt.upper()}", 0))
                    elif ref_len > alt_len:
                        alt_count = base_del
                    else:
                        alt_count = base_ins

                result = {
                    'CHROM': chrom,
                    'POS': pos,       # VCF anchor position
                    'End': end,
                    'Base_A': base_a,
                    'Base_C': base_c,
                    'Base_G': base_g,
                    'Base_T': base_t,
                    'Base_N': base_n,
                    'Base_DEL': base_del,
                    'Base_INS': base_ins,
                    'Total_Depth': total_depth,
                    'Alt_Depth': int(alt_count)
                }
                return pd.DataFrame([result])

            except Exception as ve:
                logging.warning(f"Variant processing failed for {pattern} (sample {sample}): {ve}")
                return None

        # process all variants in parallel
        with ThreadPoolExecutor(max_workers=THREADS_PER_SAMPLE) as executor:
            futures = {executor.submit(process_variant, row): idx for idx, row in df_vcf.iterrows()}
            for future in tqdm(as_completed(futures), total=len(futures), desc=f"{sample} variants"):
                try:
                    result = future.result()
                    if result is not None:
                        df_wig_combined = pd.concat([df_wig_combined, result], ignore_index=True)
                except Exception as e:
                    logging.warning(f"Thread error for sample {sample}: {e}")

        # merge only on CHROM + POS (process_variant returned one summary row per VCF POS)
        if not df_wig_combined.empty:
            merged = pd.merge(df_vcf, df_wig_combined, on=['CHROM', 'POS'], how='left')
        else:
            # if no wig info available create empty columns with 0
            merged = df_vcf.copy()
            for c in ['Base_A','Base_C','Base_G','Base_T','Base_N','Base_DEL','Base_INS','Total_Depth','Alt_Depth']:
                merged[c] = 0

        # substitute NaNs in count columns with 0 and ensure data type is integer
        count_cols = ['Base_A','Base_C','Base_G','Base_T','Base_N','Base_DEL','Base_INS','Alt_Depth']
        for c in count_cols:
            if c in merged.columns:
                merged[c] = merged[c].fillna(0).astype(int)
            else:
                merged[c] = 0

        # total_Depth: if any NaN ensure data type is int (keep 0)
        if 'Total_Depth' in merged.columns:
            merged['Total_Depth'] = merged['Total_Depth'].fillna(0).astype(int)
        else:
            merged['Total_Depth'] = 0

        # allelic burden
        def compute_allelic(row):
            td = row['Total_Depth']
            ad = row['Alt_Depth']
            if td > 0:
                return round((ad / td) * 100, 2)
            else:
                return 0.0

        merged['Allelic_Burden'] = merged.apply(compute_allelic, axis=1)

        columns_to_drop = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample_Info', 'END', 'REF_LEN', 'ALT_LEN']
        merged = merged.drop(columns=columns_to_drop, errors='ignore')

        output_dir = "Depth_Output"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{sample}_VCF_updated.csv")
        merged.to_csv(output_file, index=False)
        elapsed = round(time.time() - start_time, 2)
        logging.info(f"Sample {sample} processed and saved to {output_file} (Elapsed: {elapsed} sec)")

    except Exception as e:
        logging.error(f"Error processing {sample}: {e}")

    finally:
        # clean temporary files
        if os.path.exists(sample_temp_dir):
            for file in os.listdir(sample_temp_dir):
                try:
                    os.remove(os.path.join(sample_temp_dir, file))
                except Exception:
                    pass
            try:
                os.rmdir(sample_temp_dir)
            except Exception:
                pass
        logging.info(f"Temporary files cleaned for sample: {sample}")

def main(input_csv, n_sample_processes=4):
    df = pd.read_csv(input_csv)
    sample_list = df[['Sample_ID', 'Batch_ID', 'bedfile']].values.tolist()

    logging.info(f"Starting multiprocessing with {n_sample_processes} processes")
    with Pool(n_sample_processes) as pool:
        list(tqdm(pool.imap(process_sample, sample_list), total=len(sample_list), desc="Samples"))
    logging.info("All samples processed.")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage BamDepth_v1.py <metadata_csv> [THREADS_PER_SAMPLE]")
        logging.error("No input CSV provided. Exiting.")
    else:
        input_csv = sys.argv[1]
        # if threads provided in cmd line; update THREADS_PER_SAMPLE globally
        if len(sys.argv) > 2:
            try:
                THREADS_PER_SAMPLE = int(sys.argv[2])
            except ValueError:
                logging.warning("Invalid THREADS_PER_SAMPLE provided; using default.")
        logging.info(f"Using THREADS_PER_SAMPLE = {THREADS_PER_SAMPLE}")
        main(input_csv, n_sample_processes=4)

