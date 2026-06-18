#!/usr/bin/env python3
import csv
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import argparse


def parse_dif_file(dif_path: Path) -> list[float]:
    """Extrae los valores de ΔΔG (total energy, columna 1) del fichero Dif_ de FoldX."""
    ddg_values = []
    with open(dif_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Pdb"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    ddg_values.append(float(parts[1]))
                except ValueError:
                    pass
    return ddg_values


def prepare_foldx_and_run(mutations_csv: str, wt_pdb: str, runs: int = 5) -> None:
    df = pd.read_csv(mutations_csv)

    print(f"Repairing Wild-Type structure: {wt_pdb}...")
    subprocess.run(["foldx", "--command=RepairPDB", f"--pdb={wt_pdb}"], check=True)

    # El PDB reparado recibe el sufijo _Repair
    repaired_pdb = wt_pdb.replace(".pdb", "_Repair.pdb")
    repaired_base = Path(repaired_pdb).stem  # p. ej. physical_WT_Repair

    summary_rows = []

    for _, row in df.iterrows():
        prot_change = str(row['protein_change'])
        match = re.search(r'([A-Z])(\d+)([A-Z])', prot_change.split()[-1])

        if not match:
            continue

        wt_aa, pos, mut_aa = match.groups()
        foldx_mut = f"{wt_aa}A{pos}{mut_aa}"
        variant_name = f"{wt_aa}{pos}{mut_aa}"

        with open(f"individual_list_{pos}.txt", "w") as f:
            f.write(f"{foldx_mut};\n")

        print(f"Building mutant {variant_name} ({runs} réplicas)...")
        subprocess.run([
            "foldx",
            "--command=BuildModel",
            f"--pdb={repaired_pdb}",
            f"--mutant-file=individual_list_{pos}.txt",
            f"--numberOfRuns={runs}",
        ], check=True)

        # --- Recopilar ΔΔG de todas las réplicas ---
        # FoldX acumula filas en el mismo Dif_ si se llama varias veces desde el
        # mismo directorio; se elimina tras leerlo para evitar conteos erróneos.
        dif_file = Path(f"Dif_{repaired_base}.fxout")
        ddg_values = []
        if dif_file.exists():
            ddg_values = parse_dif_file(dif_file)
            dif_file.unlink()
        else:
            print(f"Advertencia: {dif_file} no encontrado para {variant_name}")

        # Escribir réplicas individuales
        replicas_csv = f"{variant_name}_ddg_replicas.csv"
        with open(replicas_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["replica", "ddg_kcal_mol"])
            for i, val in enumerate(ddg_values, 1):
                writer.writerow([i, val])

        # Calcular media y desviación estándar
        n = len(ddg_values)
        if n >= 2:
            mean = sum(ddg_values) / n
            sd = math.sqrt(sum((v - mean) ** 2 for v in ddg_values) / (n - 1))
        elif n == 1:
            mean = ddg_values[0]
            sd = 0.0
        else:
            mean = float("nan")
            sd = float("nan")

        summary_rows.append({
            "variant":    variant_name,
            "ddg_mean":   round(mean, 4),
            "ddg_sd":     round(sd, 4),
            "n_replicas": n,
        })

        print(f"  {variant_name}: ΔΔG = {mean:.3f} ± {sd:.3f} kcal/mol ({n} réplicas)")

        # Renombrar las N réplicas al patrón *_mutant.pdb que Nextflow espera.
        # FoldX con --numberOfRuns N genera {base}_1_0.pdb … {base}_1_{N-1}.pdb.
        for i in range(runs):
            src = Path(f"{repaired_base}_1_{i}.pdb")
            dst = Path(f"{variant_name}_rep{i}_mutant.pdb")
            if src.exists():
                shutil.move(str(src), str(dst))
                print(f"  Réplica {i} guardada como {dst.name}")
            else:
                print(f"  Advertencia: réplica {i} no encontrada ({src})")

    # Escribir resumen global de todos los mutantes
    summary_csv = "foldx_ddg_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["variant", "ddg_mean", "ddg_sd", "n_replicas"])
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nResumen ΔΔG guardado en {summary_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run FoldX Mutagenesis con réplicas.')
    parser.add_argument('-i', '--input', required=True, help='Mutations CSV')
    parser.add_argument('-p', '--pdb',   required=True, help='Cleaned WT PDB')
    parser.add_argument('-n', '--runs',  type=int, default=5,
                        help='Número de réplicas FoldX BuildModel (default: 5)')

    args = parser.parse_args()
    prepare_foldx_and_run(args.input, args.pdb, args.runs)
