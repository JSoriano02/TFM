# Predicting the Effect of DYRK1B Mutations on AZ191 Inhibitor Binding

> A reproducible **Nextflow** pipeline that combines *in silico* mutagenesis (**FoldX**) and molecular docking (**GNINA**) to predict how point mutations in the kinase **DYRK1B** affect the binding affinity of the inhibitor **AZ191**.

This repository contains the code developed for a **Master's Thesis (Trabajo de Fin de Máster, TFM)**.

---

## Table of contents

- [Overview](#overview)
- [Biological background](#biological-background)
- [Scientific approach](#scientific-approach)
- [Pipeline architecture](#pipeline-architecture)
- [Workflow diagram](#workflow-diagram)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input data](#input-data)
- [Running the pipeline](#running-the-pipeline)
- [Configuration parameters](#configuration-parameters)
- [Outputs and how to read them](#outputs-and-how-to-read-them)
- [Methodological notes](#methodological-notes)
- [Reproducibility](#reproducibility)
- [Authorship](#authorship)

---

## Overview

The goal of this project is to **predict the effect of point mutations in DYRK1B on the binding of its inhibitor AZ191**, using two complementary computational approaches that are reported separately:

1. **Structure-based *in silico* mutagenesis** with **FoldX**, which generates the structural models of each mutant and estimates the change in folding stability (ΔΔG).
2. **Molecular docking** of AZ191 onto the wild-type (WT) protein and each mutant with **GNINA**, which estimates binding affinity.

The whole workflow is orchestrated with **Nextflow (DSL2)**, providing reproducibility, modularity, caching, and the ability to resume interrupted runs (`-resume`).

A distinctive feature of this pipeline is its emphasis on **scientific rigour**: it includes a docking validation control (redocking), statistical replicates for stability estimates, fixed random seeds for reproducibility, a consistent docking-box definition shared by validation and production runs, and a protein–ligand interaction analysis to explain the results mechanistically.

---

## Biological background

**DYRK1B** (*Dual-specificity tyrosine-phosphorylation-regulated kinase 1B*) is a kinase involved in cell-cycle regulation, quiescence, and cell survival. It is a therapeutic target of interest in oncology and metabolic disease. **AZ191** is a selective small-molecule inhibitor of DYRK1B.

Mutations in the kinase domain can alter the geometry of the binding site and therefore modify inhibitor affinity — a common mechanism of drug resistance. This pipeline **quantifies that impact computationally and a priori**, distinguishing mutations that destabilise the protein fold from those that primarily disrupt inhibitor binding.

---

## Scientific approach

For each mutation the pipeline answers two independent questions:

- **Does the mutation destabilise the protein?** → measured by FoldX as ΔΔG of folding stability (kcal/mol), with 5 replicates per mutation to quantify uncertainty.
- **Does the mutation impair inhibitor binding?** → measured by GNINA as the change in docking affinity relative to the WT (kcal/mol), with the binding-pose confidence reported via the CNN score.

Crucially, these two magnitudes are **kept separate** throughout the analysis (they measure different physical properties), and each mutation is finally **classified** by which of the two it affects: stability, binding, both, or neither.

---

## Pipeline architecture

The workflow is split into independent **Nextflow processes**, defined in `main.nf` and `modules/`:

| Process | Module | Role |
|---------|--------|------|
| `FILTER_MUTATIONS` | `modules/01_filter_data.nf` | Selects the mutations of interest from the input mutation table. |
| `CLEAN_STRUCTURE` | `modules/02_clean_structure.nf` | Cleans and prepares the WT protein structure. |
| `REDOCK_VALIDATION` | `modules/00_validate_docking.nf` | **Validation control:** redocks the crystallographic AZ191 onto the WT and computes RMSD vs the native pose. |
| `FOLDX_MUTAGENESIS` | `modules/03_mutagenesis.nf` | Generates mutant structures with FoldX (5 replicates per mutation) and computes ΔΔG ± SD. |
| `PREPARE_MEEKO` | `modules/04_prepare_meeko.nf` | Prepares receptors and ligand in PDBQT format (Meeko). |
| `DOCKING_GNINA` | `modules/05_docking_gnina.nf` | Docks AZ191 onto WT and each mutant replica with GNINA, using the same box strategy as the validation. |
| `EXTRACT_AND_REPORT` | `modules/06_analysis.nf` | Aggregates affinities and CNN scores, runs the interaction analysis (ProLIF), classifies mutations, computes statistics and generates the final report. |

---

## Workflow diagram

```mermaid
flowchart TD
    A[mutations table] --> B[FILTER_MUTATIONS<br/>Select target mutations]
    C[DYRK1B-AZ191 complex PDB<br/>wild type] --> D[CLEAN_STRUCTURE<br/>Clean WT structure]

    C --> V[REDOCK_VALIDATION<br/>Redock crystal ligand + RMSD]
    D --> V
    V --> OK([Protocol validated<br/>RMSD 0.27 A])

    D --> E[FOLDX_MUTAGENESIS<br/>5 replicates per mutation<br/>delta-delta-G +/- SD]
    B --> E

    D --> F{WT + mutants}
    E --> F

    G[AZ191 ligand SDF] --> H[PREPARE_MEEKO<br/>PDBQT preparation]
    F --> H

    H --> I[DOCKING_GNINA<br/>autobox - affinity + CNN]
    I --> J[EXTRACT_AND_REPORT<br/>Aggregation - ProLIF -<br/>classification - statistics]
    D --> J

    J --> K[(Final report<br/>+ tables + plot)]
```

---

## Repository structure

```
TFM/
├── bin/                       # Python helper scripts called by the processes
│   ├── run_foldx.py           #   FoldX replicate execution + delta-delta-G aggregation
│   ├── calc_rmsd.py           #   RMSD calculation for redocking validation
│   └── generate_report.py     #   Final report, interactions, statistics
├── modules/                   # Nextflow processes (one per stage)
│   ├── 00_validate_docking.nf
│   ├── 01_filter_data.nf
│   ├── 02_clean_structure.nf
│   ├── 03_mutagenesis.nf
│   ├── 04_prepare_meeko.nf
│   ├── 05_docking_gnina.nf
│   └── 06_analysis.nf
├── raw_data/                  # Input data (NOT version-controlled)
│   ├── DYRK1B_AZ191_complex.pdb
│   └── AZ191.sdf
├── results/                   # Pipeline outputs
├── main.nf                    # Main workflow definition
├── nextflow.config            # Resources, parameters, profiles
├── .gitignore
├── README.md
└── LICENSE
```

> Script names under `bin/` are indicative; adjust them to match your repository if they differ.

---

## Requirements

| Tool | Purpose | Notes |
|------|---------|-------|
| [Nextflow](https://www.nextflow.io/) (>= 22.x, DSL2) | Workflow orchestration | Tested with 25.10.x |
| [Conda](https://docs.conda.io/) / Mamba | Environment management | `conda.enabled = true` |
| [FoldX](https://foldxsuite.crg.eu/) 4 | *In silico* mutagenesis | Academic licence; requires `rotabase.txt` |
| [GNINA](https://github.com/gnina/gnina) | Molecular docking | Requires CUDA-capable **GPU** |
| [Meeko](https://github.com/forlilab/Meeko) | Receptor/ligand preparation | Produces PDBQT |
| [ProLIF](https://prolif.readthedocs.io/) | Protein-ligand interaction fingerprints | Used in the analysis step |
| [RDKit](https://www.rdkit.org/) | Cheminformatics (RMSD, bond orders) | |
| SciPy / pandas / matplotlib / seaborn | Statistics and plotting | |
| Python >= 3.10 | Helper scripts | |

> **FoldX** requires a free academic licence and must be installed manually and available on the `PATH`. The path to `rotabase.txt` is configurable.
>
> **Graphviz** is optional; install it if you want Nextflow to render the execution DAG and HTML reports (`sudo apt install graphviz`).

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/JSoriano02/TFM.git
cd TFM

# 2. Install Nextflow if needed
curl -s https://get.nextflow.io | bash

# 3. Dependencies are resolved automatically via Conda
#    (conda.enabled = true in nextflow.config)
```

FoldX and its `rotabase.txt` must be installed separately; set the path in `nextflow.config`.

---

## Input data

The pipeline expects the following files in `raw_data/`:

| File | Description |
|------|-------------|
| `DYRK1B_AZ191_complex.pdb` | Crystallographic structure of the DYRK1B-AZ191 complex (wild type). |
| `AZ191.sdf` | Structure of the AZ191 inhibitor (ligand). |
| *mutation table* | Tabular list of mutations to evaluate. |

These files are **not included** in the repository. The crystallographic ligand is identified internally by its residue name (`QS0`); the complex also contains manganese ions (`MN`) and a phosphotyrosine residue (`PTR`), which are filtered out when the ligand is extracted.

---

## Running the pipeline

```bash
# Standard run (local executor with GPU)
nextflow run main.nf

# Resume an interrupted run without recomputing completed steps
nextflow run main.nf -resume
```

---

## Configuration parameters

Defined in `nextflow.config` and overridable from the command line (`--parameter`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `raw_dir` | `${projectDir}/raw_data` | Input data directory. |
| `outdir` | `${projectDir}/results` | Output directory. |
| `ligand_resname` | `QS0` | Residue name of the crystallographic ligand. |
| `foldx_runs` | `5` | FoldX replicates per mutation. |
| `seed` | `42` | Random seed for GNINA (reproducibility). |
| `exhaustiveness` | `16` | GNINA search exhaustiveness. |
| `num_modes` | `9` | Number of docking poses generated. |
| `rmsd_threshold` | `2.0` | Redocking validation threshold (Angstrom). |
| `stability_threshold` | `0.5` | abs(delta-delta-G) relevance threshold (kcal/mol). |
| `affinity_threshold` | `0.5` | abs(delta-affinity) relevance threshold (kcal/mol). |

**Docking box.** The search box is defined automatically via GNINA's `--autobox_ligand` (with `--autobox_add 4`) around the crystallographic ligand position. The **same box strategy is applied consistently** in both the redocking validation (`REDOCK_VALIDATION`) and the mutant docking (`DOCKING_GNINA`), so the validated protocol is exactly the one used to produce the mutant results. The legacy `box_x/y/z` and `box_size` parameters that remain in `nextflow.config` are **no longer used** by the docking processes and are kept only for reference.

> Resources (executor, CPUs, memory) and the `gpu_intensive` label (`maxForks = 1`, tuned for a 6 GB GPU) are also set in `nextflow.config`.

---

## Outputs and how to read them

After a successful run, `results/` contains:

| Output | What it tells you |
|--------|-------------------|
| `redocking_validation.txt` | RMSD of the redocked AZ191 vs the native pose. **Read this first:** the protocol passes with RMSD ≈ 0.27 Å, well below the 2 Å threshold, confirming the docking method reproduces the crystallographic pose. |
| `foldx_ddg_summary.csv` | delta-delta-G of folding stability per mutation (mean +/- SD over replicates). |
| `*_ddg_replicas.csv` | Raw per-replicate delta-delta-G values for each mutation. |
| `ddg_stability_plot.png` | Bar chart of delta-delta-G with error bars. |
| `interactions_table.csv` | ProLIF protein-ligand interactions (H-bonds, hydrophobic contacts, etc.) for WT and mutants. |
| Final report | Aggregated tables: binding affinity (GNINA Vina score **and CNN pose-confidence score**), folding stability (FoldX), mutation classification, and replicate statistics. |
| ChimeraX/PyMOL script | Visualisation of the docked poses. |

**Reading suggestion:** the WT structure is the reference row in every table. The *change* in affinity for each mutant relative to the WT is the key comparison — a positive delta-affinity means weaker binding. The **CNN score** indicates how confident the model is that a pose resembles a native binding mode (values near 1 = high confidence); it should be checked to confirm that mutant poses are reliable. Combine binding affinity with the FoldX delta-delta-G to determine whether a mutation acts on stability, on binding, or on both.

---

## Methodological notes

- **Stability vs binding are different magnitudes.** FoldX delta-delta-G quantifies folding stability; GNINA affinity quantifies ligand binding. They are reported separately and never combined into a single score.
- **Replicates and uncertainty.** FoldX is run with 5 replicates per mutation so that every delta-delta-G carries a standard deviation; one-sample statistical tests are reported against H0: delta-delta-G = 0.
- **Statistical vs biological significance.** These are treated as independent concepts: a result can be statistically significant without exceeding the biological relevance threshold (abs(delta-delta-G) > 0.5 kcal/mol), and vice versa.
- **Docking validation.** Before trusting any mutant result, the protocol is validated by redocking the native ligand and checking RMSD < 2 Å against the crystallographic pose (achieved: ≈ 0.27 Å).
- **Consistent, robust box definition.** Both the validation and the mutant docking use `--autobox_ligand` centred on the crystallographic ligand, rather than hardcoded coordinates. This guarantees that the validated protocol matches the production protocol and avoids artefacts from a free ligand centred at the origin.
- **Pose confidence.** The GNINA CNN score is reported alongside the Vina affinity so that the reliability of each docked pose can be assessed, not just its predicted energy.

---

## Reproducibility

- All stochastic steps use a fixed seed (`params.seed`); GNINA exhaustiveness and number of modes are explicit parameters.
- The full provenance (tool versions, replicates, thresholds, seed) is recorded in the final report.
- Nextflow caching allows exact resumption with `-resume`; enabling Graphviz additionally produces the execution report, timeline and DAG.

---

## Authorship

**Author:** J. Soriano ([@JSoriano02](https://github.com/JSoriano02))
**Type:** Master's Thesis (Trabajo de Fin de Máster, TFM)
**Year:** 2026

See [`LICENSE`](LICENSE) for licensing terms.