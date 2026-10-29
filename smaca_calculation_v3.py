#to run the script: python smaca_calculation_v3.py /path/to/SMACA/output/file/csv

import pandas as pd
import numpy as np
import sys
import os

def process_smaca_data(input_file, output_file):
    print(f"Reading input from: {input_file}")
    
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file was not found at:\n{input_file}")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return
    
    required_cols = ['Pi_a', 'Pi_b', 'Pi_c', 'scale_factor']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Error: The input file is missing these columns: {missing}")
        return
        
    for col in required_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    if df[required_cols].isnull().any().any():
        print("Warning: Some non-numeric values were found and converted to NaN.")

    #average of Pi_a, Pi_b, Pi_c
    #formula: =AVERAGE(B2, C2, D2)
    #AVERAGE(Pi_a, Pi_b, Pi_c)
    df['Average'] = df[['Pi_a', 'Pi_b', 'Pi_c']].mean(axis=1)

    #total copy number
    #formula: =R2*4
    #total CN: scale factor * 4
    df['Scale factor * CN'] = (df['scale_factor'] * 4).round()

    #SMN1 CN: total CN * average
    #formula: =ROUND(F2*E2)
    #Round('Scale factor * CN' * 'Average')
    df['SMN1_CN'] = (df['Scale factor * CN'] * df['Average']).round()

    #SMN2_CN: total CN - SMN1 CN
    #formula: =F2-G2
    #(Scale factor * CN) - SMN1_CN
    df['SMN2_CN'] = df['Scale factor * CN'] - df['SMN1_CN']

    #prediction (recheck/validate categorisation)
    def predict_status(smn1_count):
        if pd.isna(smn1_count):
            return "Unknown"
        
        cn = int(smn1_count)
        
        if cn <= 1:
            return "SMA"
        elif cn == 2:
            return "Normal"
        elif cn > 2:
            return "NA"
        #elif cn == 0:
            #return "SMA"
        else:
            return "NA"

    df['Prediction'] = df['SMN1_CN'].apply(predict_status)
    #reorder columns
    start_cols = ['#id', 'Pi_a', 'Pi_b', 'Pi_c']
    new_cols = ['Average', 'Scale factor * CN', 'SMN1_CN', 'SMN2_CN', 'Prediction']
    all_cols = list(df.columns)
    remaining_cols = [c for c in all_cols if c not in start_cols and c not in new_cols]
    final_order = [c for c in start_cols if c in df.columns] + new_cols + remaining_cols
    
    df_final = df[final_order]
    
    #remove footer rows in SMAca output
    if '#id' in df_final.columns:
        mask = df_final['#id'].astype(str).str.startswith('#')
        rows_dropped = mask.sum()
        if rows_dropped > 0:
            df_final = df_final[~mask]
    
    #remove columns of the 20 baseline genes in SMAca output
    if df_final.shape[1] > 21:
        df_final = df_final.iloc[:, :21]

    try:
        df_final.to_csv(output_file, index=False)
        print(f"Success! Created: {os.path.basename(output_file)}")
        print(f"Full path: {output_file}")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: You must provide a file path.")
        sys.exit(1)

    input_csv = sys.argv[1]

    if not os.path.isfile(input_csv):
        print(f"Error: Could not find the file: {input_csv}")
        sys.exit(1)
    
    input_dir = os.path.dirname(input_csv)
    input_filename = os.path.basename(input_csv)
    
    name_part, ext_part = os.path.splitext(input_filename)
    
    output_filename = f"{name_part}_calculated{ext_part}"
    output_csv = os.path.join(input_dir, output_filename)
    
    process_smaca_data(input_csv, output_csv)
