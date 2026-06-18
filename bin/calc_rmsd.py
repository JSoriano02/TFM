#!/usr/bin/env python3
"""
Computes RMSD between the best docked pose and the crystallographic reference pose.
Uses RDKit GetBestRMS to handle molecular symmetry (automorphism enumeration).
Heavy atoms only (no H), which is the standard practice in docking validation.
"""
import argparse
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


def load_template(template_sdf: str) -> Chem.Mol:
    """Loads and sanitizes the template SDF for bond-order assignment."""
    template = Chem.MolFromMolFile(template_sdf, removeHs=True)
    if template is None:
        sys.exit(f"Error: could not read template SDF: {template_sdf}")
    return template


def extract_reference_mol(ligand_ref_pdb: str, template: Chem.Mol) -> Chem.Mol:
    """
    Loads the reference ligand from the crystallographic complex PDB.
    Bond orders are assigned from the template because PDB files do not encode them.
    """
    mol_pdb = Chem.MolFromPDBFile(ligand_ref_pdb, removeHs=True, sanitize=False)
    if mol_pdb is None:
        sys.exit(f"Error: could not read reference ligand PDB: {ligand_ref_pdb}")

    try:
        mol_ref = AllChem.AssignBondOrdersFromTemplate(template, mol_pdb)
    except Exception as e:
        sys.exit(f"Error assigning bond orders to reference: {e}")

    return Chem.RemoveHs(mol_ref)


def load_best_docked_pose(docked_sdf: str, template: Chem.Mol) -> Chem.Mol:
    """
    Loads the first valid pose from a GNINA output SDF.

    GNINA writes SDFs with implicit bond orders that RDKit cannot resolve on its own
    (it may perceive incorrect valences). Strategy:
      1. Read with sanitize=False so RDKit does not reject the molecule.
      2. Keep H (removeHs=False) so AssignBondOrdersFromTemplate finds the correct atom mapping.
      3. Assign bond orders from the sanitized template.
      4. Strip H and return the clean molecule.
    Poses that fail step 3 are skipped with a warning.
    """
    suppl = Chem.SDMolSupplier(docked_sdf, sanitize=False, removeHs=False)
    for i, mol in enumerate(suppl):
        if mol is None:
            print(f"Warning: pose {i + 1} is None, skipping.")
            continue
        try:
            mol_fixed = AllChem.AssignBondOrdersFromTemplate(template, mol)
            return Chem.RemoveHs(mol_fixed)
        except Exception as e:
            print(f"Warning: pose {i + 1} failed AssignBondOrdersFromTemplate: {e}, skipping.")
            continue

    sys.exit(f"Error: no valid pose found in {docked_sdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute RMSD between the docked pose and the crystallographic reference."
    )
    parser.add_argument("--ref-pdb",   required=True, help="PDB of the ligand extracted from the complex")
    parser.add_argument("--docked",    required=True, help="SDF with GNINA output poses")
    parser.add_argument("--template",  required=True, help="Ligand template SDF (AZ191.sdf)")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="Maximum acceptable RMSD in Å (default: 2.0)")
    parser.add_argument("--output",    default="redocking_validation.txt",
                        help="Output report file")
    args = parser.parse_args()

    template   = load_template(args.template)
    mol_ref    = extract_reference_mol(args.ref_pdb, template)
    mol_docked = load_best_docked_pose(args.docked, template)

    try:
        rmsd = rdMolAlign.GetBestRMS(mol_ref, mol_docked)
    except Exception as e:
        sys.exit(f"Error computing RMSD: {e}\n"
                 "Check that the reference and the docked pose are the same ligand.")

    passed = rmsd <= args.threshold
    verdict = "PASS" if passed else "FAIL"

    lines = [
        "=== Re-docking validation (crystallographic WT) ===",
        f"RMSD best pose vs. crystal structure : {rmsd:.3f} Å",
        f"Acceptance threshold                 : {args.threshold:.1f} Å",
        f"Result                               : {verdict}",
    ]
    report = "\n".join(lines) + "\n"

    with open(args.output, "w") as f:
        f.write(report)

    print(report)

    if not passed:
        sys.exit(
            f"\n[ERROR] Docking validation FAILED: RMSD = {rmsd:.3f} Å > {args.threshold:.1f} Å.\n"
            "The docking protocol does not reproduce the crystallographic pose.\n"
            "Check the bounding box coordinates (box_x/y/z) or increase exhaustiveness."
        )


if __name__ == "__main__":
    main()
