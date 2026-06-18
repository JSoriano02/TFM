#!/usr/bin/env python3
"""
Calcula el RMSD entre la mejor pose docked y la pose cristalográfica de referencia.
Usa RDKit con GetBestRMS para manejar simetría molecular (enumeración de automorfismos).
Trabaja solo con átomos pesados (sin H), que es la práctica estándar en validación de docking.
"""
import argparse
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


def load_template(template_sdf: str) -> Chem.Mol:
    """Carga y sanitiza el SDF plantilla para asignación de órdenes de enlace."""
    template = Chem.MolFromMolFile(template_sdf, removeHs=True)
    if template is None:
        sys.exit(f"Error: no se pudo leer la plantilla SDF: {template_sdf}")
    return template


def extract_reference_mol(ligand_ref_pdb: str, template: Chem.Mol) -> Chem.Mol:
    """
    Carga el ligando de referencia desde el PDB del complejo cristalográfico.
    Asigna órdenes de enlace usando el template (los PDB no las codifican).
    """
    mol_pdb = Chem.MolFromPDBFile(ligand_ref_pdb, removeHs=True, sanitize=False)
    if mol_pdb is None:
        sys.exit(f"Error: no se pudo leer el ligando de referencia PDB: {ligand_ref_pdb}")

    try:
        mol_ref = AllChem.AssignBondOrdersFromTemplate(template, mol_pdb)
    except Exception as e:
        sys.exit(f"Error al asignar órdenes de enlace a la referencia: {e}")

    return Chem.RemoveHs(mol_ref)


def load_best_docked_pose(docked_sdf: str, template: Chem.Mol) -> Chem.Mol:
    """
    Carga la primera pose válida del SDF de salida de GNINA.

    GNINA escribe los SDF con órdenes de enlace implícitos que RDKit no puede
    resolver por sí solo (puede percibir valencias incorrectas). La estrategia:
      1. Leer con sanitize=False para que RDKit no rechace la molécula.
      2. Conservar H (removeHs=False) para que AssignBondOrdersFromTemplate
         encuentre el mapeo atómico correcto.
      3. Asignar órdenes de enlace desde el template sanitizado.
      4. Eliminar H y devolver la molécula limpia.
    Si una pose falla el paso 3, se registra y se continúa con la siguiente.
    """
    suppl = Chem.SDMolSupplier(docked_sdf, sanitize=False, removeHs=False)
    for i, mol in enumerate(suppl):
        if mol is None:
            print(f"Advertencia: pose {i + 1} es None, saltando.")
            continue
        try:
            mol_fixed = AllChem.AssignBondOrdersFromTemplate(template, mol)
            return Chem.RemoveHs(mol_fixed)
        except Exception as e:
            print(f"Advertencia: pose {i + 1} falló AssignBondOrdersFromTemplate: {e}, saltando.")
            continue

    sys.exit(f"Error: no se encontró ninguna pose válida en {docked_sdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Calcula RMSD entre pose docked y pose cristalográfica de referencia."
    )
    parser.add_argument("--ref-pdb",   required=True, help="PDB del ligando extraído del complejo")
    parser.add_argument("--docked",    required=True, help="SDF con las poses de GNINA")
    parser.add_argument("--template",  required=True, help="SDF plantilla del ligando (AZ191.sdf)")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="RMSD máximo aceptable en Å (default: 2.0)")
    parser.add_argument("--output",    default="redocking_validation.txt",
                        help="Fichero de informe de salida")
    args = parser.parse_args()

    template   = load_template(args.template)
    mol_ref    = extract_reference_mol(args.ref_pdb, template)
    mol_docked = load_best_docked_pose(args.docked, template)

    try:
        rmsd = rdMolAlign.GetBestRMS(mol_ref, mol_docked)
    except Exception as e:
        sys.exit(f"Error al calcular RMSD: {e}\n"
                 "Verifica que la referencia y la pose docked sean el mismo ligando.")

    pasa = rmsd <= args.threshold
    resultado = "PASA" if pasa else "NO PASA"

    lineas = [
        "=== Validación de redocking (WT cristalográfica) ===",
        f"RMSD mejor pose vs. cristalografía : {rmsd:.3f} Å",
        f"Umbral de aceptación               : {args.threshold:.1f} Å",
        f"Resultado                          : {resultado}",
    ]
    reporte = "\n".join(lineas) + "\n"

    with open(args.output, "w") as f:
        f.write(reporte)

    print(reporte)

    if not pasa:
        sys.exit(
            f"\n[ERROR] Validación de docking FALLIDA: RMSD = {rmsd:.3f} Å > {args.threshold:.1f} Å.\n"
            "El protocolo de docking no reproduce la pose cristalográfica.\n"
            "Revisa las coordenadas del bounding box (box_x/y/z) o aumenta exhaustiveness."
        )


if __name__ == "__main__":
    main()
