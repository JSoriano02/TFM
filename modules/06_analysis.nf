process EXTRACT_AND_REPORT {
    conda 'conda-forge::python=3.10 conda-forge::prolif conda-forge::scipy conda-forge::matplotlib'
    publishDir "${params.outdir}/06_analysis", mode: 'copy'

    input:
    path sdf_files
    path wt_pdb
    path mutant_pdbs
    path foldx_summary
    path ddg_replicas

    output:
    path "thermodynamic_report.md"
    path "visualize_interactions.cxc"
    path "interactions_table.csv"
    path "ddg_stability_plot.png"

    script:
    """
    generate_report.py \
        --sdfs ${sdf_files} \
        --wt ${wt_pdb} \
        --receptors ${mutant_pdbs} \
        --foldx-summary ${foldx_summary} \
        --replicas ${ddg_replicas} \
        --stability-threshold ${params.ddg_stability_threshold} \
        --affinity-threshold ${params.delta_affinity_threshold} \
        --foldx-runs ${params.foldx_runs} \
        --seed ${params.seed} \
        --exhaustiveness ${params.exhaustiveness} \
        --rmsd-threshold ${params.rmsd_threshold}
    """
}