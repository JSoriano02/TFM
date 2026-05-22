#!/usr/bin/env python3
import pandas as pd
import argparse
import subprocess
import re
import sys
import shutil
from pathlib import Path

def prepare_foldx_and_run(mutations_csv: str, wt_pdb: str) -> None:
    df = pd.read_csv(mutations_csv)
    
    print(f"Repairing Wild-Type structure: {wt_pdb}...")
    subprocess.run(["foldx", "--command=RepairPDB", f"--pdb={wt_pdb}"], check=True)
    
    # The repaired PDB is named by appending _Repair
    repaired_pdb = wt_pdb.replace(".pdb", "_Repair.pdb")
    
    for _, row in df.iterrows():
        prot_change = str(row['protein_change'])
        match = re.search(r'([A-Z])(\d+)([A-Z])', prot_change.split()[-1])
        
        if match:
            wt_aa, pos, mut_aa = match.groups()
            foldx_mut = f"{wt_aa}A{pos}{mut_aa}"
            
            # Create a clean identifier for the final file (e.g., E112K)
            variant_name = f"{wt_aa}{pos}{mut_aa}"
            
            with open(f"individual_list_{pos}.txt", "w") as f:
                f.write(f"{foldx_mut};\n")
            
            print(f"Building mutant {variant_name}...")
            subprocess.run([
                "foldx", 
                "--command=BuildModel", 
                f"--pdb={repaired_pdb}", 
                f"--mutant-file=individual_list_{pos}.txt"
            ], check=True)
            
            # FoldX always generates the mutated structure as {repaired_pdb}_1.pdb
            expected_out = repaired_pdb.replace(".pdb", "_1.pdb")
            final_name = f"{variant_name}_mutant.pdb"
            
            # Rename the file immediately to prevent it from being overwritten
            # by the next iteration of the loop
            if Path(expected_out).exists():
                shutil.move(expected_out, final_name)
                print(f"Successfully saved {final_name}")
            else:
                print(f"Warning: FoldX output not found for {variant_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run FoldX Mutagenesis.')
    parser.add_argument('-i', '--input', required=True, help='Mutations CSV')
    parser.add_argument('-p', '--pdb', required=True, help='Cleaned WT PDB')
    
    args = parser.parse_args()
    prepare_foldx_and_run(args.input, args.pdb)