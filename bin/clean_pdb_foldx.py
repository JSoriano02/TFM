#!/usr/bin/env python3
"""
clean_pdb_foldx.py
Strict PDB cleaner designed specifically to prevent FoldX 4 parsing crashes.
"""

import argparse
import sys

def clean_for_foldx(input_pdb: str, output_pdb: str) -> None:
    # Standard amino acids accepted by FoldX
    valid_res = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
                 
    # Map common modified/phosphorylated residues to their standard counterparts
    mod_res = {"SEP": "SER", "TPO": "THR", "PTR": "TYR", "MSE": "MET", "CSO": "CYS"}
    
    first_chain = None
    
    try:
        with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
            for line in fin:
                if line.startswith("ATOM  "):
                    # 1. Keep only the primary alternate conformation (' ' or 'A' or '1')
                    alt_loc = line[16]
                    if alt_loc not in [' ', 'A', '1']:
                        continue
                        
                    # 2. Isolate the first chain only to avoid multimer clashes
                    chain = line[21]
                    if first_chain is None:
                        first_chain = chain
                    if chain != first_chain:
                        continue
                        
                    # 3. Fix modified residues (e.g., Phosphotyrosine to Tyrosine)
                    resname = line[17:20].strip()
                    if resname in mod_res:
                        new_res = mod_res[resname]
                        line = line[:17] + new_res + line[20:]
                        resname = new_res
                        
                    # 4. Skip any remaining non-standard molecules
                    if resname not in valid_res:
                        continue
                        
                    # 5. Erase the alt_loc flag to prevent FoldX from getting confused
                    line = line[:16] + ' ' + line[17:]
                    
                    # 6. CRITICAL FIX: Ensure strict 80-character line width
                    line = line.rstrip('\r\n')
                    line = line.ljust(80) + '\n'
                    
                    fout.write(line)
                    
            # Conclude the PDB file correctly
            fout.write("TER\nEND\n")
            
    except IOError as error:
        print(f"Error processing PDB: {error}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Strict PDB cleaner for FoldX compatibility")
    parser.add_argument('-i', '--input', required=True, help="Input raw PDB file")
    parser.add_argument('-o', '--output', required=True, help="Output cleaned PDB file")
    args = parser.parse_args()
    
    clean_for_foldx(args.input, args.output)