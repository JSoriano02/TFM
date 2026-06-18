#!/usr/bin/env python3

import argparse
import re
import subprocess
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# External tool version detection
# ---------------------------------------------------------------------------

def get_tool_versions() -> Dict[str, str]:
    """Queries the actual versions of FoldX and GNINA at runtime."""
    versions: Dict[str, str] = {}

    # FoldX
    try:
        result = subprocess.run(
            ["foldx"], capture_output=True, text=True, timeout=10
        )
        output = result.stdout + result.stderr
        match = re.search(r'FoldX\s+(\d[\w.]*)', output)
        versions["foldx"] = f"FoldX {match.group(1)}" if match else "FoldX 4 (exact version unavailable)"
    except Exception:
        versions["foldx"] = "FoldX (not found in PATH)"

    # GNINA
    try:
        result = subprocess.run(
            ["gnina", "--version"], capture_output=True, text=True, timeout=10
        )
        output = (result.stdout + result.stderr).strip()
        versions["gnina"] = output if output else "version unavailable"
    except Exception:
        versions["gnina"] = "GNINA (not found in PATH)"

    return versions


# ---------------------------------------------------------------------------
# GNINA score parsing
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
# ChimeraX visualisation script
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
# Mutation classification: stability vs. binding
# ---------------------------------------------------------------------------

def _normalize_variant(name: str) -> str:
    """Strips the '_mutant' suffix to cross-reference GNINA names with FoldX names."""
    return name.replace("_mutant", "")


def _base_mutation(name: str) -> str:
    """Strips replica and mutant suffixes to extract the base mutation name.
    Examples: E112K_rep0_mutant → E112K
              DYRK1B_AZ191_complex_clean → DYRK1B_AZ191_complex_clean  (WT, unchanged)
    """
    name = name.replace("_mutant", "")
    return re.sub(r"_rep\d+$", "", name)


def aggregate_docking_results(results: Dict[str, dict]) -> Dict[str, dict]:
    """
    Groups per-replica docking scores by base mutation name.
    WT has a single docking run; mutants have one run per FoldX replica.
    """
    from collections import defaultdict
    groups: dict = defaultdict(list)
    for key, data in results.items():
        groups[_base_mutation(key)].append(data["vina"])

    aggregated = {}
    for base, vina_list in groups.items():
        n = len(vina_list)
        mean_vina = sum(vina_list) / n
        sd_vina = (sum((v - mean_vina)**2 for v in vina_list) / (n - 1))**0.5 if n > 1 else 0.0
        first_key = next(k for k in results if _base_mutation(k) == base)
        aggregated[base] = {
            "vina":      round(mean_vina, 3),
            "vina_sd":   round(sd_vina, 3),
            "n_docked":  n,
            "cnn":       results[first_key]["cnn"],
            "file_path": results[first_key]["file_path"],
        }
    return aggregated


def classify_mutations(
    results: Dict[str, dict],
    foldx_df: Optional[pd.DataFrame],
    stability_threshold: float,
    affinity_threshold: float,
) -> Optional[pd.DataFrame]:
    """
    Builds a classification table by cross-referencing FoldX (stability) and
    GNINA (binding affinity) data. The two quantities are always kept separate.

    Output columns: variant, ddg_mean, ddg_sd, vina, delta_affinity, classification
    """
    if foldx_df is None:
        return None

    # Identify the WT: the SDF that does not match any FoldX mutant
    foldx_variants = set(foldx_df["variant"].astype(str))
    wt_key = next(
        (k for k in results if _base_mutation(k) not in foldx_variants),
        None,
    )
    wt_vina = results[wt_key]["vina"] if wt_key else 0.0

    rows = []
    for result_key, data in results.items():
        norm = _base_mutation(result_key)
        is_wt = (result_key == wt_key)

        if is_wt:
            ddg_mean, ddg_sd = 0.0, 0.0
            delta_aff = 0.0
            classification = "WT (reference)"
        else:
            foldx_row = foldx_df[foldx_df["variant"].astype(str) == norm]
            if foldx_row.empty:
                print(f"Warning: '{norm}' not found in foldx_summary. ΔΔG will be NaN.")
                ddg_mean, ddg_sd = float("nan"), float("nan")
            else:
                ddg_mean = float(foldx_row["ddg_mean"].iloc[0])
                ddg_sd   = float(foldx_row["ddg_sd"].iloc[0])

            delta_aff = data["vina"] - wt_vina

            affects_stability = abs(ddg_mean)  > stability_threshold
            affects_binding   = abs(delta_aff) > affinity_threshold

            if affects_stability and affects_binding:
                classification = "Both"
            elif affects_stability:
                classification = "Stability"
            elif affects_binding:
                classification = "Binding"
            else:
                classification = "Neither"

        rows.append({
            "variant":        result_key,
            "ddg_stab_mean":  ddg_mean,
            "ddg_stab_sd":    ddg_sd,
            "vina_kcal_mol":  data["vina"],
            "delta_affinity": delta_aff,
            "classification": classification,
        })

    return pd.DataFrame(rows).set_index("variant")


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
    Generates the Markdown report with six sections:
      0. Traceability (versions, date, key parameters)
      1. Binding affinity (GNINA)
      2. Folding stability (FoldX)
      3. Per-mutant classification
      4. Statistics over FoldX replicas
      5. ProLIF interactions (if available)
    """
    output_report = Path("thermodynamic_report.md")
    sorted_results = sorted(results.items(), key=lambda item: item[1]['vina'])

    try:
        with output_report.open("w") as md:

            md.write("# Analysis report: DYRK1B / AZ191\n\n")

            # ------------------------------------------------------------------
            # Section 0: Traceability
            # ------------------------------------------------------------------
            if run_params:
                tool_v = run_params.get("tool_versions", {})
                md.write("## Traceability\n\n")
                md.write("| Parameter | Value |\n")
                md.write("|-----------|-------|\n")
                md.write(f"| Run date                  | {date.today()} |\n")
                md.write(f"| FoldX                     | {tool_v.get('foldx', 'N/A')} |\n")
                md.write(f"| GNINA                     | {tool_v.get('gnina', 'N/A')} |\n")
                md.write(f"| FoldX replicas            | {run_params.get('foldx_runs', 'N/A')} |\n")
                md.write(f"| GNINA seed                | {run_params.get('seed', 'N/A')} |\n")
                md.write(f"| GNINA exhaustiveness      | {run_params.get('exhaustiveness', 'N/A')} |\n")
                md.write(f"| Validation RMSD threshold | {run_params.get('rmsd_threshold', 'N/A')} Å |\n")
                md.write(f"| ΔΔG stability threshold   | {run_params.get('stability_threshold', 'N/A')} kcal/mol |\n")
                md.write(f"| Δ affinity threshold      | {run_params.get('affinity_threshold', 'N/A')} kcal/mol |\n")
                md.write("\n")

            md.write("> **Methodological note:** FoldX measures ΔΔG of *folding stability* "
                     "of the protein; GNINA measures *binding affinity* of the ligand. "
                     "These are distinct quantities and are **always reported separately**.\n\n")

            # ------------------------------------------------------------------
            # Section 1: Binding affinity (GNINA)
            # ------------------------------------------------------------------
            md.write("## 1. Binding affinity — GNINA (kcal/mol)\n\n")
            md.write("More negative Vina score indicates higher predicted affinity.\n\n")
            md.write("| Mutation | Mean Vina (kcal/mol) | ± SD | n |\n")
            md.write("|----------|---------------------|------|---|\n")
            for name, data in sorted_results:
                sd_str = f"±{data['vina_sd']:.3f}" if data['n_docked'] > 1 else "—"
                md.write(f"| {name} | {data['vina']:.2f} | {sd_str} | {data['n_docked']} |\n")

            md.write("\n*Vina score*: more negative = higher predicted affinity. "
                     "SD across docked replicas; WT has a single docking run.\n\n")

            # ------------------------------------------------------------------
            # Section 2: Folding stability (FoldX)
            # ------------------------------------------------------------------
            md.write("## 2. Folding stability — FoldX (kcal/mol)\n\n")
            if classification_df is not None:
                md.write("ΔΔG > 0 indicates destabilisation. "
                         f"Relevance threshold: |ΔΔG| > {stability_threshold} kcal/mol.\n\n")
                md.write("| Variant | Mean ΔΔG (kcal/mol) | ± SD |\n")
                md.write("|---------|---------------------|------|\n")
                for variant, row in classification_df.iterrows():
                    mean = row["ddg_stab_mean"]
                    sd   = row["ddg_stab_sd"]
                    if variant == next(
                        (k for k in results if row["classification"] == "WT (reference)"),
                        None,
                    ):
                        md.write(f"| {variant} | 0.00 (reference) | — |\n")
                    else:
                        md.write(f"| {variant} | {mean:+.3f} | ±{sd:.3f} |\n")
            else:
                md.write("*FoldX data not available.*\n")

            # ------------------------------------------------------------------
            # Section 3: Classification
            # ------------------------------------------------------------------
            md.write("\n## 3. Mutation classification\n\n")
            if classification_df is not None:
                md.write(
                    f"Thresholds: |ΔΔG_FoldX| > {stability_threshold} kcal/mol "
                    f"(stability) and |Δaffinity| > {affinity_threshold} kcal/mol (binding).\n\n"
                )
                md.write("| Variant | ΔΔG stab. (kcal/mol) | Δ affinity (kcal/mol) | Category |\n")
                md.write("|---------|----------------------|-----------------------|---------|\n")
                for variant, row in classification_df.iterrows():
                    mean      = row["ddg_stab_mean"]
                    delta_aff = row["delta_affinity"]
                    cat       = row["classification"]
                    if cat == "WT (reference)":
                        md.write(f"| {variant} | 0.00 | 0.00 | {cat} |\n")
                    else:
                        md.write(f"| {variant} | {mean:+.3f} | {delta_aff:+.3f} | {cat} |\n")
            else:
                md.write("*FoldX data not available for classification.*\n")

            # ------------------------------------------------------------------
            # Section 4: FoldX replica statistics
            # ------------------------------------------------------------------
            if stats_df is not None and not stats_df.empty:
                md.write("\n## 4. FoldX replica statistics\n\n")
                md.write(
                    f"Tests against H₀: ΔΔG = 0 (WT as reference). "
                    f"Biological relevance threshold: |ΔΔG| > {stability_threshold} kcal/mol.\n"
                    "> **Note:** statistical significance and biological relevance are "
                    "independent concepts. A result can be statistically significant "
                    "without being biologically relevant, and vice versa.\n\n"
                )
                md.write("| Variant | n | Mean ΔΔG ± SD | Test | p-value | Sig. (p<0.05) | Biol. relevant |\n")
                md.write("|---------|---|--------------|------|---------|--------------|---------------|\n")
                for variant, row in stats_df.iterrows():
                    p_str  = f"{row['p_value']:.4f}" if not pd.isna(row["p_value"]) else "N/A"
                    sig    = "✓" if row["significant_p005"]      else "–"
                    rel    = "✓" if row["biologically_relevant"] else "–"
                    md.write(
                        f"| {variant} | {int(row['n_replicas'])} | "
                        f"{row['ddg_mean']:+.3f} ± {row['ddg_sd']:.3f} | "
                        f"{row['test']} | {p_str} | {sig} | {rel} |\n"
                    )
                md.write("\n*See `ddg_stability_plot.png` for error bars.*\n")

            # ------------------------------------------------------------------
            # Section 5: ProLIF interactions
            # ------------------------------------------------------------------
            if interactions_df is not None and not interactions_df.empty:
                md.write("\n## 5. Protein-ligand interactions (ProLIF)\n\n")
                md.write("Comparative table WT vs. mutants. ✓ = interaction present; – = absent.\n\n")

                cols = interactions_df.columns.tolist()
                header = "| Variant | " + " | ".join(
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
# FoldX replica statistics
# ---------------------------------------------------------------------------

def compute_statistics(
    replica_csvs: List[str],
    stability_threshold: float,
) -> Optional[pd.DataFrame]:
    """
    For each FoldX replica CSV ({variant}_ddg_replicas.csv):
      - Shapiro-Wilk normality test.
      - One-sample t-test against 0 (if normal) or Wilcoxon signed-rank (if not normal).
      - Distinguishes statistical significance (p < 0.05) from biological relevance
        (|mean ΔΔG| > stability_threshold).
    """
    try:
        from scipy import stats as spstats
    except ImportError:
        print("Warning: scipy not available. Skipping statistical analysis.")
        return None

    rows = []
    for csv_path_str in replica_csvs:
        csv_path = Path(csv_path_str)
        variant = csv_path.stem.replace("_ddg_replicas", "")

        try:
            df_rep = pd.read_csv(csv_path)
            values = df_rep["ddg_kcal_mol"].dropna().values
        except Exception as e:
            print(f"Warning: could not read {csv_path}: {e}")
            continue

        n = len(values)
        if n == 0:
            continue

        mean = float(values.mean())
        sd   = float(values.std(ddof=1)) if n > 1 else 0.0

        # Normality test requires at least 3 data points
        if n >= 3:
            _, p_shapiro = spstats.shapiro(values)
            normal = bool(p_shapiro > 0.05)
        else:
            normal = True

        # Test against H₀: μ = 0
        if n >= 2 and normal:
            _, p_test = spstats.ttest_1samp(values, 0.0)
            test_name = "t-test (1 sample)"
        elif n >= 2:
            try:
                _, p_test = spstats.wilcoxon(values)
                test_name = "Wilcoxon"
            except ValueError:
                p_test = float("nan")
                test_name = "Wilcoxon (N/A)"
        else:
            p_test = float("nan")
            test_name = "insufficient N"

        significant = (not pd.isna(p_test)) and (p_test < 0.05)
        relevant    = abs(mean) > stability_threshold

        rows.append({
            "variant":               variant,
            "n_replicas":            n,
            "ddg_mean":              round(mean, 4),
            "ddg_sd":                round(sd, 4),
            "test":                  test_name,
            "p_value":               round(p_test, 4) if not pd.isna(p_test) else float("nan"),
            "significant_p005":      significant,
            "biologically_relevant": relevant,
        })

    if not rows:
        return None

    return pd.DataFrame(rows).set_index("variant")


def generate_ddg_plot(
    stats_df: pd.DataFrame,
    stability_threshold: float,
) -> None:
    """
    Bar chart of mean ΔΔG ± SD per mutant.
    Red = statistically significant (p < 0.05), blue = not significant.
    Reference lines at 0 and ±stability_threshold.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available. Skipping plot generation.")
        return

    variants = stats_df.index.tolist()
    means    = stats_df["ddg_mean"].values
    sds      = stats_df["ddg_sd"].values
    sigs     = stats_df["significant_p005"].values

    colors = ["#d62728" if s else "#1f77b4" for s in sigs]

    fig, ax = plt.subplots(figsize=(max(6, len(variants) * 1.8), 5))

    bars = ax.bar(variants, means, yerr=sds, color=colors,
                  capsize=6, alpha=0.82, edgecolor="black", linewidth=0.8)

    ax.axhline(y=0, color="black", linewidth=1.2, linestyle="-", label="WT (reference)")
    ax.axhline(y= stability_threshold, color="gray", linewidth=1.0, linestyle="--",
               label=f"Threshold ±{stability_threshold} kcal/mol")
    ax.axhline(y=-stability_threshold, color="gray", linewidth=1.0, linestyle="--")

    # Asterisk above significant bars
    for bar, sd_val, sig in zip(bars, sds, sigs):
        if sig:
            y_pos = bar.get_height() + sd_val + abs(max(means) - min(means)) * 0.03
            ax.text(bar.get_x() + bar.get_width() / 2, y_pos, "*",
                    ha="center", va="bottom", fontsize=16, color="black")

    ax.set_xlabel("Variant", fontsize=12)
    ax.set_ylabel("ΔΔG folding (kcal/mol)", fontsize=12)
    ax.set_title("Effect of mutations on DYRK1B stability (FoldX)", fontsize=13)
    ax.legend(fontsize=10)
    ax.tick_params(axis="x", labelsize=11)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d62728", label="Significant (p < 0.05)"),
        Patch(facecolor="#1f77b4", label="Not significant"),
    ]
    ax.legend(handles=legend_elements + ax.get_legend_handles_labels()[0][:2], fontsize=10)

    plt.tight_layout()
    plt.savefig("ddg_stability_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Plot saved to ddg_stability_plot.png")


# ---------------------------------------------------------------------------
# ProLIF interaction analysis
# ---------------------------------------------------------------------------

def analyze_prolif_interactions(
    sdf_files: List[str],
    wt_pdb: str,
    mutant_pdbs: List[str],
) -> Optional[pd.DataFrame]:
    """
    Computes protein-ligand interaction fingerprints with ProLIF for each variant.
    Returns a DataFrame with rows = variants, columns = (residue, interaction_type).
    Returns None if ProLIF is unavailable or fails for all variants.
    """
    try:
        import MDAnalysis as mda
        import prolif as plf
    except ImportError as e:
        print(f"Warning: could not import ProLIF/MDAnalysis ({e}). "
              "Skipping interaction analysis.")
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
            print(f"Warning: no receptor PDB found for '{variant}'. "
                  "Skipping interactions for this variant.")
            continue

        try:
            u_prot = mda.Universe(receptor_pdb)
            protein_mol = plf.Molecule.from_mda(u_prot.select_atoms("protein"))

            # sanitize=False to tolerate the implicit valences that GNINA writes in the SDF
            suppl = _Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
            rdkit_lig = next((m for m in suppl if m is not None), None)
            if rdkit_lig is None:
                raise ValueError(f"no valid poses in {sdf_path}")
            rdkit_lig.UpdatePropertyCache(strict=False)
            _Chem.FastFindRings(rdkit_lig)
            ligand_mol = plf.Molecule.from_rdkit(rdkit_lig)

            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand_mol], protein_mol)

            df = fp.to_dataframe()
            if not df.empty:
                frames[variant] = df.iloc[0]

        except Exception as exc:
            print(f"Warning: ProLIF failed for '{variant}': {exc}")
            continue

    # Always write the CSV even with no results: Nextflow declares it a mandatory output
    result = pd.DataFrame(frames).T if frames else pd.DataFrame()
    result.index.name = "variant"
    result.to_csv("interactions_table.csv")

    if frames:
        print(f"Interactions table saved to interactions_table.csv "
              f"({len(result)} variants, {len(result.columns)} interactions)")
    else:
        print("Warning: ProLIF produced no results. interactions_table.csv written empty.")

    return result if not result.empty else None


# ---------------------------------------------------------------------------
# Entry point
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

    foldx_df: Optional[pd.DataFrame] = None
    if foldx_summary_path and Path(foldx_summary_path).exists():
        try:
            foldx_df = pd.read_csv(foldx_summary_path)
        except Exception as e:
            print(f"Warning: could not read {foldx_summary_path}: {e}")
    else:
        print("Warning: --foldx-summary not provided or file not found. "
              "Stability section will be unavailable.")

    tool_versions      = get_tool_versions()
    aggregated_results = aggregate_docking_results(results)
    interactions_df    = analyze_prolif_interactions(sdf_files, wt_pdb, mutant_pdbs)
    classification_df  = classify_mutations(
        aggregated_results, foldx_df, stability_threshold, affinity_threshold
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
        aggregated_results, interactions_df, classification_df, stats_df,
        stability_threshold, affinity_threshold, run_params,
    )

    print("Analysis complete: thermodynamic report, ChimeraX script, "
          "interactions table, classification, statistics and traceability generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract GNINA scores, analyse ProLIF interactions and classify mutations."
    )
    parser.add_argument('--sdfs',       nargs='+', required=True,
                        help="List of docked SDF files")
    parser.add_argument('--wt',         required=True,
                        help="Cleaned Wild-Type PDB structure")
    parser.add_argument('--receptors',  nargs='*', default=[],
                        help="PDB files of FoldX-generated mutants")
    parser.add_argument('--foldx-summary', default=None,
                        help="CSV with mean ΔΔG and SD from FoldX (foldx_ddg_summary.csv)")
    parser.add_argument('--replicas',   nargs='*', default=[],
                        help="CSVs with individual FoldX replicas ({variant}_ddg_replicas.csv)")
    parser.add_argument('--stability-threshold', type=float, default=0.5,
                        help="|ΔΔG_FoldX| threshold for classifying a stability change (kcal/mol)")
    parser.add_argument('--affinity-threshold',  type=float, default=0.5,
                        help="|Δaffinity_GNINA| threshold for classifying a binding change (kcal/mol)")
    # Traceability parameters
    parser.add_argument('--foldx-runs',      type=int,   default=5,
                        help="Number of FoldX replicas used in the run")
    parser.add_argument('--seed',            type=int,   default=42,
                        help="Seed used for GNINA")
    parser.add_argument('--exhaustiveness',  type=int,   default=16,
                        help="Exhaustiveness used for GNINA")
    parser.add_argument('--rmsd-threshold',  type=float, default=2.0,
                        help="RMSD threshold from docking validation (Å)")

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
