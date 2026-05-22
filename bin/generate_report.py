#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Tuple, Dict, List

def parse_sdf_scores(sdf_path: Path) -> Tuple[float, float]:
    """
    Extracts Vina affinity and CNN scores from the first pose of a Gnina output SDF.
    Uses robust regex to handle inconsistent whitespace and line endings.
    """
    vina_score = 0.0
    cnn_score = 0.0
    
    try:
        content = sdf_path.read_text(encoding='utf-8')
        
        # \s* maneja cualquier tipo de salto de línea o espacio
        # ([-\d.]+) captura explícitamente números positivos, negativos y decimales
        vina_match = re.search(r'<minimizedAffinity>\s*\n([-\d.]+)', content)
        if vina_match:
            vina_score = float(vina_match.group(1))
            
        cnn_match = re.search(r'<CNNscore>\s*\n([-\d.]+)', content)
        if cnn_match:
            cnn_score = float(cnn_match.group(1))
            
    except IOError as error:
        print(f"Error reading file {sdf_path}: {error}")
            
    return vina_score, cnn_score


def generate_chimerax_script(results: Dict[str, dict], wt_pdb_name: str) -> None:
    """
    Generates a ChimeraX command script (.cxc) to visualize the interactions.
    
    Args:
        results (Dict): Dictionary containing docking scores and file paths.
        wt_pdb_name (str): Filename of the wild-type PDB structure.
    """
    output_script = Path("visualize_interactions.cxc")
    
    try:
        with output_script.open("w") as cxc:
            cxc.write("# ChimeraX script for AZ191 - DYRK1B visualization\n")
            cxc.write("set bgColor white\n\n")
            
            # 1. Load the Wild-Type Receptor (Assigned as Model #1)
            cxc.write(f"open ../02_cleaned/{wt_pdb_name}\n")
            
            # 2. Load all docked poses of the mutants (Models #2 onwards)
            for data in results.values():
                cxc.write(f"open {data['file_path']}\n")
                
            # 3. Clean up the default view
            cxc.write("\nhide atoms\n")
            
            # 4. Display the receptor as ribbons and apply color
            cxc.write("show cartoons #1\n")
            cxc.write("color #1 slate\n")
            
            # 5. Display the docked ligands as sticks
            cxc.write("\n# Format docked ligands\n")
            cxc.write("show atoms #2-\n")
            cxc.write("color #2- byhetero\n")
            cxc.write("color #2- green target c\n")
            
            # 6. Calculate and display hydrogen bonds visually
            cxc.write("\n# Compute H-bonds between receptor (#1) and ligands (#2-)\n")
            cxc.write("hbonds #1 restrict #2- reveal true color black lineThickness 3\n")
            
            # 7. Focus the camera on the binding site
            cxc.write("view #2-\n")
            
    except IOError as error:
        print(f"Error writing ChimeraX script: {error}")


def generate_thermodynamic_report(results: Dict[str, dict]) -> None:
    """
    Generates a Markdown report justifying the thermodynamic results.
    
    Args:
        results (Dict): Dictionary containing docking scores and file paths.
    """
    output_report = Path("thermodynamic_report.md")
    
    # Sort results by Vina Affinity (most negative / favorable first)
    sorted_results = sorted(results.items(), key=lambda item: item[1]['vina'])
    
    try:
        with output_report.open("w") as md:
            md.write("# Thermodynamic Justification of AZ191 Binding\n\n")
            
            md.write("## Overview of Empirical Scoring\n")
            md.write("Binding affinity ($\\Delta G_{bind}$) is approximated by the docking score. ")
            md.write("According to the Gibbs free energy equation ($\\Delta G = \\Delta H - T\\Delta S$), ")
            md.write("a more negative Vina score indicates an energetically favorable complex, primarily driven ")
            md.write("by enthalpic contributions ($\\Delta H$) such as hydrogen bonding and van der Waals ")
            md.write("interactions within the kinase hinge region.\n\n")
            
            md.write("## Results\n")
            md.write("| Structure | Vina Affinity (kcal/mol) | CNN Score (0-1) |\n")
            md.write("|-----------|--------------------------|-----------------|\n")
            
            for name, data in sorted_results:
                md.write(f"| {name} | {data['vina']:.2f} | {data['cnn']:.3f} |\n")
                
            md.write("\n## Discussion\n")
            md.write("1. **Vina Affinity**: Represents the standard empirical free energy. Mutations ")
            md.write("altering steric bulk or electrostatic properties in the binding pocket will shift this score.\n")
            md.write("2. **CNN Score**: Indicates the probability (0 to 1) that the generated pose is ")
            md.write("highly accurate (RMSD < 2Å to native). A drop in CNN score in mutant variants ")
            md.write("relative to the wild-type indicates that the somatic mutation may structurally ")
            md.write("disrupt the primary binding mode of AZ191.\n")
            
    except IOError as error:
        print(f"Error writing thermodynamic report: {error}")


def main(sdf_files: List[str], wt_pdb: str) -> None:
    """
    Main execution workflow.
    """
    results = {}
    
    for sdf_str in sdf_files:
        sdf_path = Path(sdf_str)
        # Assuming filename format like 'E112K_docked.sdf'
        variant_name = sdf_path.stem.replace('_docked', '') 
        
        vina_val, cnn_val = parse_sdf_scores(sdf_path)
        
        results[variant_name] = {
            'vina': vina_val,
            'cnn': cnn_val,
            'file_path': sdf_str
        }
        
    generate_chimerax_script(results, wt_pdb)
    generate_thermodynamic_report(results)
    
    print("Analysis complete: ChimeraX script and thermodynamic report generated successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Gnina scores and prepare visualization.")
    parser.add_argument('--sdfs', nargs='+', required=True, help="List of docked SDF files")
    parser.add_argument('--wt', required=True, help="Filename of the Wild-Type PDB structure")
    
    args = parser.parse_args()
    main(args.sdfs, args.wt)