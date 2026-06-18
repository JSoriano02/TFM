#!/usr/bin/env python3

import argparse
import re
import subprocess
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Trazabilidad: versiones de herramientas externas (Prioridad 6)
# ---------------------------------------------------------------------------

def get_tool_versions() -> Dict[str, str]:
    """Consulta las versiones reales de FoldX y GNINA en tiempo de ejecución."""
    versions: Dict[str, str] = {}

    # FoldX
    try:
        result = subprocess.run(
            ["foldx"], capture_output=True, text=True, timeout=10
        )
        output = result.stdout + result.stderr
        match = re.search(r'FoldX\s+(\d[\w.]*)', output)
        versions["foldx"] = f"FoldX {match.group(1)}" if match else "FoldX 4 (versión exacta no disponible)"
    except Exception:
        versions["foldx"] = "FoldX (no encontrado en PATH)"

    # GNINA
    try:
        result = subprocess.run(
            ["gnina", "--version"], capture_output=True, text=True, timeout=10
        )
        output = (result.stdout + result.stderr).strip()
        versions["gnina"] = output if output else "versión no disponible"
    except Exception:
        versions["gnina"] = "GNINA (no encontrado en PATH)"

    return versions


# ---------------------------------------------------------------------------
# Parseo de scores GNINA
# ---------------------------------------------------------------------------

def parse_sdf_scores(sdf_path: Path) -> Tuple[float, float]:
    """
    Extracts Vina affinity and CNN scores from the first pose of a Gnina output SDF.
    Uses robust regex to handle inconsistent whitespace and line endings.
    """
    vina_score = 0.0
    cnn_score = 0.0

    try:
        content = sdf_path.read_text(encoding='utf-8')

        vina_match = re.search(r'<minimizedAffinity>\s*\n([-\d.]+)', content)
        if vina_match:
            vina_score = float(vina_match.group(1))

        cnn_match = re.search(r'<CNNscore>\s*\n([-\d.]+)', content)
        if cnn_match:
            cnn_score = float(cnn_match.group(1))

    except IOError as error:
        print(f"Error reading file {sdf_path}: {error}")

    return vina_score, cnn_score


# ---------------------------------------------------------------------------
# Script ChimeraX (sin cambios respecto a la versión anterior)
# ---------------------------------------------------------------------------

def generate_chimerax_script(results: Dict[str, dict], wt_pdb_name: str) -> None:
    """Generates a ChimeraX command script (.cxc) to visualize the interactions."""
    output_script = Path("visualize_interactions.cxc")

    try:
        with output_script.open("w") as cxc:
            cxc.write("# ChimeraX script for AZ191 - DYRK1B visualization\n")
            cxc.write("set bgColor white\n\n")
            cxc.write(f"open ../02_cleaned/{wt_pdb_name}\n")
            for data in results.values():
                cxc.write(f"open {data['file_path']}\n")
            cxc.write("\nhide atoms\n")
            cxc.write("show cartoons #1\n")
            cxc.write("color #1 slate\n")
            cxc.write("\n# Format docked ligands\n")
            cxc.write("show atoms #2-\n")
            cxc.write("color #2- byhetero\n")
            cxc.write("color #2- green target c\n")
            cxc.write("\n# Compute H-bonds between receptor (#1) and ligands (#2-)\n")
            cxc.write("hbonds #1 restrict #2- reveal true color black lineThickness 3\n")
            cxc.write("view #2-\n")
    except IOError as error:
        print(f"Error writing ChimeraX script: {error}")


# ---------------------------------------------------------------------------
# Clasificación estabilidad vs unión (Prioridad 4)
# ---------------------------------------------------------------------------

def _normalize_variant(name: str) -> str:
    """Elimina el sufijo '_mutant' para cruzar nombres GNINA con nombres FoldX."""
    return name.replace("_mutant", "")


def _base_mutation(name: str) -> str:
    """Extrae la mutación base descartando sufijos de réplica y de mutant.
    Ejemplos: E112K_rep0_mutant → E112K
              DYRK1B_AZ191_complex_clean → DYRK1B_AZ191_complex_clean  (WT, sin cambio)
    """
    name = name.replace("_mutant", "")
    return re.sub(r"_rep\d+$", "", name)


def classify_mutations(
    results: Dict[str, dict],
    foldx_df: Optional[pd.DataFrame],
    stability_threshold: float,
    affinity_threshold: float,
) -> Optional[pd.DataFrame]:
    """
    Construye una tabla de clasificación cruzando datos de FoldX (estabilidad) y
    GNINA (afinidad de unión). Magnitudes tratadas siempre por separado.

    Columnas de salida: variante, ddg_mean, ddg_sd, vina, delta_affinity, clasificacion
    """
    if foldx_df is None:
        return None

    # Identificar la WT: el SDF que no corresponde a ningún mutante FoldX
    foldx_variants = set(foldx_df["variant"].astype(str))
    wt_key = next(
        (k for k in results if _base_mutation(k) not in foldx_variants),
        None,
    )
    wt_vina = results[wt_key]["vina"] if wt_key else 0.0

    rows = []
    for result_key, data in results.items():
        norm = _base_mutation(result_key)   # E112K_rep0_mutant → E112K
        is_wt = (result_key == wt_key)

        if is_wt:
            ddg_mean, ddg_sd = 0.0, 0.0
            delta_aff = 0.0
            clasificacion = "WT (referencia)"
        else:
            foldx_row = foldx_df[foldx_df["variant"].astype(str) == norm]
            if foldx_row.empty:
                print(f"Advertencia: '{norm}' no encontrado en foldx_summary. "
                      "Se usará ΔΔG = NaN.")
                ddg_mean, ddg_sd = float("nan"), float("nan")
            else:
                ddg_mean = float(foldx_row["ddg_mean"].iloc[0])
                ddg_sd   = float(foldx_row["ddg_sd"].iloc[0])

            delta_aff = data["vina"] - wt_vina

            afecta_estab  = abs(ddg_mean)   > stability_threshold
            afecta_union  = abs(delta_aff)  > affinity_threshold

            if afecta_estab and afecta_union:
                clasificacion = "Ambas"
            elif afecta_estab:
                clasificacion = "Estabilidad"
            elif afecta_union:
                clasificacion = "Unión"
            else:
                clasificacion = "Ninguna"

        rows.append({
            "variante":       result_key,
            "ΔΔG_estab_mean": ddg_mean,
            "ΔΔG_estab_sd":   ddg_sd,
            "vina_kcal_mol":  data["vina"],
            "Δ_afinidad":     delta_aff,
            "clasificacion":  clasificacion,
        })

    return pd.DataFrame(rows).set_index("variante")


# ---------------------------------------------------------------------------
# Informe Markdown
# ---------------------------------------------------------------------------

def generate_thermodynamic_report(
    results: Dict[str, dict],
    interactions_df: Optional[pd.DataFrame],
    classification_df: Optional[pd.DataFrame],
    stats_df: Optional[pd.DataFrame],
    stability_threshold: float,
    affinity_threshold: float,
    run_params: Optional[Dict] = None,
) -> None:
    """
    Genera el informe Markdown con seis secciones:
      0. Trazabilidad (versiones, fecha, parámetros clave)
      1. Afinidad de unión (GNINA)
      2. Estabilidad de plegamiento (FoldX)
      3. Clasificación por mutante
      4. Estadística sobre réplicas FoldX
      5. Interacciones ProLIF (si están disponibles)
    """
    output_report = Path("thermodynamic_report.md")
    sorted_results = sorted(results.items(), key=lambda item: item[1]['vina'])

    try:
        with output_report.open("w") as md:

            md.write("# Informe de análisis: DYRK1B / AZ191\n\n")

            # ------------------------------------------------------------------
            # Sección 0: Trazabilidad
            # ------------------------------------------------------------------
            if run_params:
                tool_v = run_params.get("tool_versions", {})
                md.write("## Trazabilidad\n\n")
                md.write("| Parámetro | Valor |\n")
                md.write("|-----------|-------|\n")
                md.write(f"| Fecha de ejecución       | {date.today()} |\n")
                md.write(f"| FoldX                    | {tool_v.get('foldx', 'N/A')} |\n")
                md.write(f"| GNINA                    | {tool_v.get('gnina', 'N/A')} |\n")
                md.write(f"| Réplicas FoldX           | {run_params.get('foldx_runs', 'N/A')} |\n")
                md.write(f"| Semilla GNINA            | {run_params.get('seed', 'N/A')} |\n")
                md.write(f"| Exhaustiveness GNINA     | {run_params.get('exhaustiveness', 'N/A')} |\n")
                md.write(f"| Umbral RMSD validación   | {run_params.get('rmsd_threshold', 'N/A')} Å |\n")
                md.write(f"| Umbral ΔΔG estabilidad   | {run_params.get('stability_threshold', 'N/A')} kcal/mol |\n")
                md.write(f"| Umbral Δ afinidad        | {run_params.get('affinity_threshold', 'N/A')} kcal/mol |\n")
                md.write("\n")

            md.write("> **Nota metodológica:** FoldX mide ΔΔG de *estabilidad de plegamiento* "
                     "de la proteína; GNINA mide *afinidad de unión* del ligando. "
                     "Son magnitudes distintas y se reportan **siempre por separado**.\n\n")

            # ------------------------------------------------------------------
            # Sección 1: Afinidad de unión (GNINA)
            # ------------------------------------------------------------------
            md.write("## 1. Afinidad de unión — GNINA (kcal/mol)\n\n")
            md.write("Cuanto más negativo el valor Vina, mayor afinidad predicha.\n\n")
            md.write("| Estructura | Vina (kcal/mol) | CNN Score |\n")
            md.write("|------------|-----------------|----------|\n")
            for name, data in sorted_results:
                md.write(f"| {name} | {data['vina']:.2f} | {data['cnn']:.3f} |\n")

            md.write("\n*CNN Score*: probabilidad de que la pose tenga RMSD < 2 Å respecto "
                     "a la pose nativa. Caída en mutantes indica disrupción del modo de unión.\n\n")

            # ------------------------------------------------------------------
            # Sección 2: Estabilidad de plegamiento (FoldX)
            # ------------------------------------------------------------------
            md.write("## 2. Estabilidad de plegamiento — FoldX (kcal/mol)\n\n")
            if classification_df is not None:
                md.write("ΔΔG > 0 indica desestabilización. "
                         f"Umbral de relevancia: |ΔΔG| > {stability_threshold} kcal/mol.\n\n")
                md.write("| Variante | ΔΔG medio (kcal/mol) | ± SD |\n")
                md.write("|----------|---------------------|------|\n")
                for variant, row in classification_df.iterrows():
                    mean = row["ΔΔG_estab_mean"]
                    sd   = row["ΔΔG_estab_sd"]
                    if variant == next(
                        (k for k in results if row["clasificacion"] == "WT (referencia)"),
                        None,
                    ):
                        md.write(f"| {variant} | 0.00 (referencia) | — |\n")
                    else:
                        md.write(f"| {variant} | {mean:+.3f} | ±{sd:.3f} |\n")
            else:
                md.write("*Datos de FoldX no disponibles.*\n")

            # ------------------------------------------------------------------
            # Sección 3: Clasificación
            # ------------------------------------------------------------------
            md.write("\n## 3. Clasificación de mutaciones\n\n")
            if classification_df is not None:
                md.write(
                    f"Umbrales: |ΔΔG_FoldX| > {stability_threshold} kcal/mol "
                    f"(estabilidad) y |Δafinidad| > {affinity_threshold} kcal/mol (unión).\n\n"
                )
                md.write("| Variante | ΔΔG estab. (kcal/mol) | Δ afinidad (kcal/mol) | Categoría |\n")
                md.write("|----------|-----------------------|-----------------------|----------|\n")
                for variant, row in classification_df.iterrows():
                    mean      = row["ΔΔG_estab_mean"]
                    delta_aff = row["Δ_afinidad"]
                    cat       = row["clasificacion"]
                    if cat == "WT (referencia)":
                        md.write(f"| {variant} | 0.00 | 0.00 | {cat} |\n")
                    else:
                        md.write(f"| {variant} | {mean:+.3f} | {delta_aff:+.3f} | {cat} |\n")
            else:
                md.write("*Datos de FoldX no disponibles para clasificación.*\n")

            # ------------------------------------------------------------------
            # Sección 4: Estadística (réplicas FoldX)
            # ------------------------------------------------------------------
            if stats_df is not None and not stats_df.empty:
                md.write("\n## 4. Estadística sobre réplicas FoldX\n\n")
                md.write(
                    f"Tests realizados contra H₀: ΔΔG = 0 (WT como referencia). "
                    f"Umbral de relevancia biológica: |ΔΔG| > {stability_threshold} kcal/mol.\n"
                    "> **Nota:** significancia estadística y relevancia biológica son "
                    "conceptos independientes. Un resultado puede ser estadísticamente "
                    "significativo sin ser biológicamente relevante, y viceversa.\n\n"
                )
                md.write("| Variante | n | ΔΔG medio ± SD | Test | p-valor | Sig. (p<0.05) | Relevante biol. |\n")
                md.write("|----------|---|---------------|------|---------|--------------|----------------|\n")
                for variant, row in stats_df.iterrows():
                    p_str  = f"{row['p_valor']:.4f}" if not pd.isna(row["p_valor"]) else "N/A"
                    sig    = "✓" if row["significativo_p005"]      else "–"
                    rel    = "✓" if row["relevante_biologicamente"] else "–"
                    md.write(
                        f"| {variant} | {int(row['n_replicas'])} | "
                        f"{row['ddg_mean']:+.3f} ± {row['ddg_sd']:.3f} | "
                        f"{row['test']} | {p_str} | {sig} | {rel} |\n"
                    )
                md.write("\n*El gráfico `ddg_stability_plot.png` muestra las barras de error.*\n")

            # ------------------------------------------------------------------
            # Sección 5: Interacciones ProLIF
            # ------------------------------------------------------------------
            if interactions_df is not None and not interactions_df.empty:
                md.write("\n## 5. Interacciones proteína-ligando (ProLIF)\n\n")
                md.write("Tabla comparativa WT vs mutantes. ✓ = interacción presente; – = ausente.\n\n")

                cols = interactions_df.columns.tolist()
                header = "| Variante | " + " | ".join(
                    f"{res} {itype}" for res, itype in cols
                ) + " |"
                separator = "|---|" + "|".join(["---"] * len(cols)) + "|"
                md.write(header + "\n")
                md.write(separator + "\n")
                for variant, row in interactions_df.iterrows():
                    cells = " | ".join("✓" if v else "–" for v in row)
                    md.write(f"| {variant} | {cells} |\n")

    except IOError as error:
        print(f"Error writing thermodynamic report: {error}")


# ---------------------------------------------------------------------------
# Estadística sobre réplicas FoldX (Prioridad 5)
# ---------------------------------------------------------------------------

def compute_statistics(
    replica_csvs: List[str],
    stability_threshold: float,
) -> Optional[pd.DataFrame]:
    """
    Para cada CSV de réplicas FoldX ({variant}_ddg_replicas.csv):
      - Test de normalidad Shapiro-Wilk.
      - t-test de una muestra contra 0 (si normal) o Wilcoxon signed-rank (si no normal).
      - Distingue significancia estadística (p < 0.05) de relevancia biológica
        (|ΔΔG_medio| > stability_threshold).
    """
    try:
        from scipy import stats as spstats
    except ImportError:
        print("Advertencia: scipy no disponible. Saltando análisis estadístico.")
        return None

    rows = []
    for csv_path_str in replica_csvs:
        csv_path = Path(csv_path_str)
        variant = csv_path.stem.replace("_ddg_replicas", "")

        try:
            df_rep = pd.read_csv(csv_path)
            values = df_rep["ddg_kcal_mol"].dropna().values
        except Exception as e:
            print(f"Advertencia: no se pudo leer {csv_path}: {e}")
            continue

        n = len(values)
        if n == 0:
            continue

        mean = float(values.mean())
        sd   = float(values.std(ddof=1)) if n > 1 else 0.0

        # Test de normalidad (mínimo 3 puntos)
        if n >= 3:
            _, p_shapiro = spstats.shapiro(values)
            normal = bool(p_shapiro > 0.05)
        else:
            normal = True

        # Test contra H₀: μ = 0
        if n >= 2 and normal:
            _, p_test = spstats.ttest_1samp(values, 0.0)
            test_name = "t-test (1 muestra)"
        elif n >= 2:
            try:
                _, p_test = spstats.wilcoxon(values)
                test_name = "Wilcoxon"
            except ValueError:
                p_test = float("nan")
                test_name = "Wilcoxon (N/A)"
        else:
            p_test = float("nan")
            test_name = "N insuficiente"

        significant = (not pd.isna(p_test)) and (p_test < 0.05)
        relevant    = abs(mean) > stability_threshold

        rows.append({
            "variante":               variant,
            "n_replicas":             n,
            "ddg_mean":               round(mean, 4),
            "ddg_sd":                 round(sd, 4),
            "test":                   test_name,
            "p_valor":                round(p_test, 4) if not pd.isna(p_test) else float("nan"),
            "significativo_p005":     significant,
            "relevante_biologicamente": relevant,
        })

    if not rows:
        return None

    return pd.DataFrame(rows).set_index("variante")


def generate_ddg_plot(
    stats_df: pd.DataFrame,
    stability_threshold: float,
) -> None:
    """
    Gráfico de barras ΔΔG medio ± SD para cada mutante.
    Rojo = estadísticamente significativo (p < 0.05), azul = no significativo.
    Líneas de referencia en 0 y en ±stability_threshold.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Advertencia: matplotlib no disponible. Saltando generación del gráfico.")
        return

    variants = stats_df.index.tolist()
    means    = stats_df["ddg_mean"].values
    sds      = stats_df["ddg_sd"].values
    sigs     = stats_df["significativo_p005"].values

    colors = ["#d62728" if s else "#1f77b4" for s in sigs]

    fig, ax = plt.subplots(figsize=(max(6, len(variants) * 1.8), 5))

    bars = ax.bar(variants, means, yerr=sds, color=colors,
                  capsize=6, alpha=0.82, edgecolor="black", linewidth=0.8)

    ax.axhline(y=0, color="black", linewidth=1.2, linestyle="-", label="WT (referencia)")
    ax.axhline(y= stability_threshold, color="gray", linewidth=1.0, linestyle="--",
               label=f"Umbral ±{stability_threshold} kcal/mol")
    ax.axhline(y=-stability_threshold, color="gray", linewidth=1.0, linestyle="--")

    # Asterisco sobre barras significativas
    for bar, sd_val, sig in zip(bars, sds, sigs):
        if sig:
            y_pos = bar.get_height() + sd_val + abs(max(means) - min(means)) * 0.03
            ax.text(bar.get_x() + bar.get_width() / 2, y_pos, "*",
                    ha="center", va="bottom", fontsize=16, color="black")

    ax.set_xlabel("Variante", fontsize=12)
    ax.set_ylabel("ΔΔG plegamiento (kcal/mol)", fontsize=12)
    ax.set_title("Efecto de las mutaciones sobre la estabilidad de DYRK1B (FoldX)", fontsize=13)
    ax.legend(fontsize=10)
    ax.tick_params(axis="x", labelsize=11)

    # Leyenda de colores
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d62728", label="Significativo (p < 0.05)"),
        Patch(facecolor="#1f77b4", label="No significativo"),
    ]
    ax.legend(handles=legend_elements + ax.get_legend_handles_labels()[0][:2], fontsize=10)

    plt.tight_layout()
    plt.savefig("ddg_stability_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Gráfico guardado en ddg_stability_plot.png")


# ---------------------------------------------------------------------------
# Análisis de interacciones ProLIF (sin cambios respecto a Prioridad 3)
# ---------------------------------------------------------------------------

def analyze_prolif_interactions(
    sdf_files: List[str],
    wt_pdb: str,
    mutant_pdbs: List[str],
) -> Optional[pd.DataFrame]:
    """
    Calcula fingerprints de interacción proteína-ligando con ProLIF para cada variante.
    Devuelve DataFrame con filas = variantes, columnas = (residuo, tipo_interacción).
    Si ProLIF no está disponible o falla, devuelve None.
    """
    try:
        import MDAnalysis as mda
        import prolif as plf
    except ImportError as e:
        print(f"Advertencia: no se pudo importar ProLIF/MDAnalysis ({e}). "
              "Saltando análisis de interacciones.")
        return None

    from rdkit import Chem as _Chem

    receptor_map: Dict[str, str] = {}
    for pdb in [wt_pdb] + list(mutant_pdbs):
        receptor_map[Path(pdb).stem] = pdb

    frames: Dict[str, pd.Series] = {}

    for sdf_str in sdf_files:
        sdf_path = Path(sdf_str)
        variant = sdf_path.stem.replace("_docked", "")
        receptor_pdb = receptor_map.get(variant)

        if receptor_pdb is None:
            print(f"Advertencia: no se encontró PDB receptor para '{variant}'. "
                  "Saltando interacciones de esta variante.")
            continue

        try:
            # Proteína: MDAnalysis sobre PDB (sin cambios)
            u_prot = mda.Universe(receptor_pdb)
            protein_mol = plf.Molecule.from_mda(u_prot.select_atoms("protein"))

            # Ligando: RDKit con sanitize=False para tolerar las valencias
            # implícitas que escribe GNINA en el SDF.
            suppl = _Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
            rdkit_lig = next((m for m in suppl if m is not None), None)
            if rdkit_lig is None:
                raise ValueError(f"sin poses válidas en {sdf_path}")
            rdkit_lig.UpdatePropertyCache(strict=False)
            _Chem.FastFindRings(rdkit_lig)
            ligand_mol = plf.Molecule.from_rdkit(rdkit_lig)

            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand_mol], protein_mol)

            df = fp.to_dataframe()
            if not df.empty:
                frames[variant] = df.iloc[0]

        except Exception as exc:
            print(f"Advertencia: ProLIF falló para '{variant}': {exc}")
            continue

    # Siempre crear el CSV aunque no haya resultados: Nextflow lo declara output obligatorio
    result = pd.DataFrame(frames).T if frames else pd.DataFrame()
    result.index.name = "variante"
    result.to_csv("interactions_table.csv")

    if frames:
        print(f"Tabla de interacciones guardada en interactions_table.csv "
              f"({len(result)} variantes, {len(result.columns)} interacciones)")
    else:
        print("Advertencia: ProLIF no produjo resultados. interactions_table.csv guardado vacío.")

    return result if not result.empty else None


# ---------------------------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------------------------

def main(
    sdf_files: List[str],
    wt_pdb: str,
    mutant_pdbs: List[str],
    foldx_summary_path: Optional[str],
    replica_csvs: List[str],
    stability_threshold: float,
    affinity_threshold: float,
    foldx_runs: int,
    seed: int,
    exhaustiveness: int,
    rmsd_threshold: float,
) -> None:

    results = {}
    for sdf_str in sdf_files:
        sdf_path = Path(sdf_str)
        variant_name = sdf_path.stem.replace('_docked', '')
        vina_val, cnn_val = parse_sdf_scores(sdf_path)
        results[variant_name] = {
            'vina': vina_val,
            'cnn': cnn_val,
            'file_path': sdf_str,
        }

    # Cargar resumen FoldX
    foldx_df: Optional[pd.DataFrame] = None
    if foldx_summary_path and Path(foldx_summary_path).exists():
        try:
            foldx_df = pd.read_csv(foldx_summary_path)
        except Exception as e:
            print(f"Advertencia: no se pudo leer {foldx_summary_path}: {e}")
    else:
        print("Advertencia: --foldx-summary no proporcionado o fichero no encontrado. "
              "Sección de estabilidad no disponible.")

    tool_versions     = get_tool_versions()
    interactions_df   = analyze_prolif_interactions(sdf_files, wt_pdb, mutant_pdbs)
    classification_df = classify_mutations(
        results, foldx_df, stability_threshold, affinity_threshold
    )
    stats_df = compute_statistics(replica_csvs, stability_threshold)

    if stats_df is not None and not stats_df.empty:
        generate_ddg_plot(stats_df, stability_threshold)

    run_params = {
        "tool_versions":      tool_versions,
        "foldx_runs":         foldx_runs,
        "seed":               seed,
        "exhaustiveness":     exhaustiveness,
        "rmsd_threshold":     rmsd_threshold,
        "stability_threshold": stability_threshold,
        "affinity_threshold":  affinity_threshold,
    }

    generate_chimerax_script(results, wt_pdb)
    generate_thermodynamic_report(
        results, interactions_df, classification_df, stats_df,
        stability_threshold, affinity_threshold, run_params,
    )

    print("Análisis completo: informe termodinámico, script ChimeraX, "
          "tabla de interacciones, clasificación, estadística y trazabilidad generados.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extrae scores GNINA, analiza interacciones ProLIF y clasifica mutaciones."
    )
    parser.add_argument('--sdfs',       nargs='+', required=True,
                        help="Lista de ficheros SDF docked")
    parser.add_argument('--wt',         required=True,
                        help="PDB de la estructura Wild-Type limpia")
    parser.add_argument('--receptors',  nargs='*', default=[],
                        help="PDB de los mutantes generados por FoldX")
    parser.add_argument('--foldx-summary', default=None,
                        help="CSV con ΔΔG medios y SD de FoldX (foldx_ddg_summary.csv)")
    parser.add_argument('--replicas',   nargs='*', default=[],
                        help="CSVs con réplicas individuales FoldX ({variant}_ddg_replicas.csv)")
    parser.add_argument('--stability-threshold', type=float, default=0.5,
                        help="Umbral |ΔΔG_FoldX| para clasificar cambio en estabilidad (kcal/mol)")
    parser.add_argument('--affinity-threshold',  type=float, default=0.5,
                        help="Umbral |Δafinidad_GNINA| para clasificar cambio en unión (kcal/mol)")
    # Parámetros de trazabilidad (Prioridad 6)
    parser.add_argument('--foldx-runs',      type=int,   default=5,
                        help="Número de réplicas FoldX usadas en la ejecución")
    parser.add_argument('--seed',            type=int,   default=42,
                        help="Semilla usada en GNINA")
    parser.add_argument('--exhaustiveness',  type=int,   default=16,
                        help="Exhaustiveness usada en GNINA")
    parser.add_argument('--rmsd-threshold',  type=float, default=2.0,
                        help="Umbral RMSD de la validación de docking (Å)")

    args = parser.parse_args()
    main(
        args.sdfs,
        args.wt,
        args.receptors,
        args.foldx_summary,
        args.replicas,
        args.stability_threshold,
        args.affinity_threshold,
        args.foldx_runs,
        args.seed,
        args.exhaustiveness,
        args.rmsd_threshold,
    )
