process EXTRACT_AND_REPORT {
    conda 'conda-forge::python=3.10'
    publishDir "${params.outdir}/06_analysis", mode: 'copy'

    input:
    path sdf_files
    path wt_pdb

    output:
    path "thermodynamic_report.md"
    // Updated from .pml to .cxc to match the ChimeraX script generation
    path "visualize_interactions.cxc"

    script:
    """
    generate_report.py --sdfs ${sdf_files} --wt ${wt_pdb}
    """
}