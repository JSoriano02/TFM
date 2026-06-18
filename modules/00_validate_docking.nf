process REDOCK_VALIDATION {
    label 'gpu_intensive'
    conda 'conda-forge::meeko conda-forge::pdb-tools'
    publishDir "${params.outdir}/00_validation", mode: 'copy'

    input:
    path wt_pdb        // cleaned WT structure (output of CLEAN_STRUCTURE)
    path complex_pdb   // original crystallographic complex (to extract the reference pose)
    path ligand_sdf    // AZ191.sdf (ligand structure)

    output:
    path "redocking_validation.txt", emit: report

    script:
    """
    # 1. Prepare WT receptor in PDBQT format
    pdb_element ${wt_pdb} > fixed_wt.pdb
    mk_prepare_receptor.py -i fixed_wt.pdb -o wt_receptor -p -a

    # 2. Prepare ligand in PDBQT format
    mk_prepare_ligand.py -i ${ligand_sdf} -o ligand.pdbqt

    # 3. Extract the crystallographic ligand from the complex (non-water HETATM).
    #    Must be done BEFORE docking: used as autobox_ligand to center the search box
    #    on the real binding site, and as the RMSD reference pose.
    grep -E "^HETATM" ${complex_pdb} | grep " ${params.ligand_resname} " > crystal_ligand.pdb
    echo "END" >> crystal_ligand.pdb

    # 4. Re-dock AZ191 onto the WT.
    #    --autobox_ligand centers the search box on the crystallographic ligand coordinates,
    #    avoiding the issue of the free AZ191.sdf being centered at the origin (≈0,0,0).
    gnina -r wt_receptor.pdbqt -l ligand.pdbqt \
        --autobox_ligand crystal_ligand.pdb \
        --autobox_add 4 \
        --cnn_scoring rescore \
        --seed ${params.seed} \
        --exhaustiveness ${params.exhaustiveness} \
        --num_modes ${params.num_modes} \
        --out redocked_wt.sdf

    # 5. Compute RMSD and emit report (exits with error if RMSD > threshold)
    calc_rmsd.py \
        --ref-pdb  crystal_ligand.pdb \
        --docked   redocked_wt.sdf \
        --template ${ligand_sdf} \
        --threshold ${params.rmsd_threshold} \
        --output   redocking_validation.txt
    """
}
