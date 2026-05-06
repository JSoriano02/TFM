#!/usr/bin/env python3
import pandas as pd
import argparse
import re

def filter_gdc_mutations(input_tsv, output_csv, top_n=3):
    # 1. Load the raw GDC TSV file
    df = pd.read_csv(input_tsv, sep='\t')
    
    # 2. Filter strictly for 'Missense' mutations using the exact GDC column name
    df_missense = df[df['consequence'].str.contains('Missense', case=False, na=False)].copy()
    
    # 3. Extract the amino acid position from the 'protein_change' column
    def extract_position(change_str):
        parts = str(change_str).split()
        if len(parts) > 1:
            # Search for the first sequence of digits in the second word 
            match = re.search(r'\d+', parts[1])
            return int(match.group()) if match else -1
        return -1
        
    df_missense['AA_Position'] = df_missense['protein_change'].apply(extract_position)
    
    # 4. Filter by the DYRK1B kinase domain boundaries (residues 78 to 442)
    df_domain = df_missense[(df_missense['AA_Position'] >= 78) & (df_missense['AA_Position'] <= 442)].copy()
    
    # 5. Extract the number of affected patients using the exact GDC column name
    df_domain['Patient_Count'] = df_domain['num_ssm_affected_cases'].fillna(0).astype(int)
    
    # 6. Sort by patient count in descending order
    df_sorted = df_domain.sort_values(by='Patient_Count', ascending=False)
    
    # 7. Keep the top N most frequent mutations for the docking pipeline
    df_final = df_sorted.head(top_n)
    
    # 8. Save a clean CSV file with the essential columns 
    columns_to_keep = [
        'ssm_id', 'protein_change', 'consequence', 
        'num_ssm_affected_cases', 'sift_impact', 'polyphen_impact'
    ]
    
    # Ensure we only select columns that actually exist in the dataframe
    final_cols = [c for c in columns_to_keep if c in df_final.columns]
    
    df_final[final_cols].to_csv(output_csv, index=False)
    
    print(f"Filtering complete: Isolated the top {len(df_final)} missense mutations within the kinase domain.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter GDC mutations for docking.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file from GDC')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    parser.add_argument('-n', '--top', type=int, default=3, help='Maximum number of mutations to keep')
    args = parser.parse_args()
    
    filter_gdc_mutations(args.input, args.output, args.top)