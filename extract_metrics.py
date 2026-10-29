#to run the script: python extract_metrics.py /path/to/sample/metrics/directory/

import os
import glob
import pandas as pd
import io
import sys

def get_picard_metric_row(filepath, category_col=None, category_val=None):
    if not os.path.exists(filepath):
        return None
    try:
        clean_lines = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('## HISTOGRAM'):
                    break
                if not line.startswith('#') and line.strip():
                    clean_lines.append(line)
        if not clean_lines:
            return None

        df = pd.read_csv(io.StringIO(''.join(clean_lines)), sep='\t')
        
        if df.empty:
            return None

        if category_col and category_val:
            if category_col in df.columns:
                subset = df[df[category_col] == category_val]
                if not subset.empty:
                    return subset.iloc[0].to_dict()
            return df.iloc[0].to_dict()
        
        return df.iloc[0].to_dict()

    except Exception as e:
        print(f"Warning: Could not parse {filepath}: {e}")
        return None

def main():
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = "."

    if not os.path.isdir(target_dir):
        print(f"Error: Directory not found: {target_dir}")
        return

    print(f"Scanning directory: {target_dir}")

    search_pattern = os.path.join(target_dir, "*_alignment_metrics.txt")
    align_files = glob.glob(search_pattern)
    
    if not align_files:
        print(f"No *_alignment_metrics.txt files found in {target_dir}")
        return

    print(f"Found {len(align_files)} samples. Processing...")
    
    collected_data = []

    for al_file_path in align_files:
        base_filename = os.path.basename(al_file_path)
        sample_name = base_filename.replace("_alignment_metrics.txt", "")
        
        markdup_file = os.path.join(target_dir, f"{sample_name}_metrics.txt")
        insert_file = os.path.join(target_dir, f"{sample_name}_insert_size_metrics.txt")
        wgs_file = os.path.join(target_dir, f"{sample_name}_wgs_metrics.txt")
        
        row = {'Sample': sample_name}

        al_data = get_picard_metric_row(al_file_path, category_col='CATEGORY', category_val='PAIR')
        if al_data:
            row['TOTAL_READS'] = al_data.get('TOTAL_READS', '')
            row['PF_READS_ALIGNED'] = al_data.get('PF_READS_ALIGNED', '')
            row['MEDIAN_READ_LENGTH'] = al_data.get('MEDIAN_READ_LENGTH', '')
            row['STRAND_BALANCE'] = al_data.get('STRAND_BALANCE', '')

        md_data = get_picard_metric_row(markdup_file)
        if md_data:
            row['PERCENT_DUPLICATION'] = md_data.get('PERCENT_DUPLICATION', '')
            row['ESTIMATED_LIBRARY_SIZE'] = md_data.get('ESTIMATED_LIBRARY_SIZE', '')
            row['UNMAPPED_READS'] = md_data.get('UNMAPPED_READS', '')

        is_data = get_picard_metric_row(insert_file)
        if is_data:
            row['MEDIAN_INSERT_SIZE'] = is_data.get('MEDIAN_INSERT_SIZE', '')

        wgs_data = get_picard_metric_row(wgs_file)
        if wgs_data:
            row['MEAN_COVERAGE'] = wgs_data.get('MEAN_COVERAGE', '')
            row['SD_COVERAGE'] = wgs_data.get('SD_COVERAGE', '')
            row['MEDIAN_COVERAGE'] = wgs_data.get('MEDIAN_COVERAGE', '')
            
            pct_fields = ['PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_15X', 'PCT_20X', 
                          'PCT_25X', 'PCT_30X', 'PCT_40X', 'PCT_50X']
            for field in pct_fields:
                row[field] = wgs_data.get(field, '')

        collected_data.append(row)

    df_final = pd.DataFrame(collected_data)
    
    cols_order = [
        'Sample',
        'TOTAL_READS', 'PF_READS_ALIGNED', 'MEDIAN_READ_LENGTH', 'STRAND_BALANCE',
        'PERCENT_DUPLICATION', 'ESTIMATED_LIBRARY_SIZE', 'UNMAPPED_READS',
        'MEDIAN_INSERT_SIZE',
        'MEAN_COVERAGE', 'SD_COVERAGE', 'MEDIAN_COVERAGE', 
        'PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_15X', 'PCT_20X', 
        'PCT_25X', 'PCT_30X', 'PCT_40X', 'PCT_50X'
    ]
    
    final_cols = [c for c in cols_order if c in df_final.columns]
    df_final = df_final[final_cols]

    output_filename = 'combined_metrics_summary.tsv'
    df_final.to_csv(output_filename, sep='\t', index=False)
    print(f"Success! Aggregated metrics saved to: {os.path.abspath(output_filename)}")

if __name__ == "__main__":
    main()
