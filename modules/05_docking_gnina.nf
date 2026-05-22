process DOCKING_GNINA {
    label 'gpu_intensive'
    publishDir "${params.outdir}/05_docking", mode: 'copy'

    input:
    path receptor_pdbqt
    path ligand_pdbqt

    output:
    path "${receptor_pdbqt.baseName}_docked.sdf", emit: docked_sdf

    script:
    """
    # Gnina is automatically found in the bin/ directory.
    # CNN scoring explicitly set to 'rescore' to balance speed and accuracy on the RTX 2060
    gnina -r ${receptor_pdbqt} -l ${ligand_pdbqt} \
        --center_x ${params.box_x} --center_y ${params.box_y} --center_z ${params.box_z} \
        --size_x ${params.box_size} --size_y ${params.box_size} --size_z ${params.box_size} \
        --cnn_scoring rescore \
        --exhaustiveness 8 \
        --out ${receptor_pdbqt.baseName}_docked.sdf
    """
}