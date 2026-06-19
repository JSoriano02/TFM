process DOCKING_GNINA {
    label 'gpu_intensive'
    publishDir "${params.outdir}/05_docking", mode: 'copy'

    input:
    path receptor_pdbqt
    path ligand_pdbqt
    path complex_pdb    // crystallographic complex; used only to extract the reference binding site

    output:
    path "${receptor_pdbqt.baseName}_docked.sdf", emit: docked_sdf

    script:
    """
    # Extract the crystallographic ligand to define the search box.
    # Identical strategy to REDOCK_VALIDATION: centres the box on the real binding site
    # coordinates rather than relying on manually measured params.box_x/y/z.
    grep -E "^HETATM" ${complex_pdb} | grep " ${params.ligand_resname} " > crystal_ligand.pdb
    echo "END" >> crystal_ligand.pdb

    gnina -r ${receptor_pdbqt} -l ${ligand_pdbqt} \\
        --autobox_ligand crystal_ligand.pdb \\
        --autobox_add 4 \\
        --cnn_scoring rescore \\
        --seed ${params.seed} \\
        --exhaustiveness ${params.exhaustiveness} \\
        --num_modes ${params.num_modes} \\
        --out ${receptor_pdbqt.baseName}_docked.sdf
    """
}