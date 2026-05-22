process PREPARE_MEEKO {
    // Added pdb-tools to the conda environment to fix FoldX PDB formatting
    conda 'conda-forge::meeko conda-forge::pdb-tools'
    publishDir "${params.outdir}/04_prepared", mode: 'copy'

    input:
    path receptor_pdb
    path ligand_sdf

    output:
    path "${receptor_pdb.baseName}.pdbqt", emit: receptor_pdbqt
    path "${ligand_sdf.baseName}.pdbqt", emit: ligand_pdbqt

    script:
    """
    # 1. FoldX strips the element column (77-78) which crashes Meeko/RDKit.
    # We use pdb_element to reconstruct the missing element symbols based on atom names.
    pdb_element ${receptor_pdb} > fixed_receptor.pdb

    # 2. Prepare receptor allowing incomplete surface residues (-a)
    # Feed the fixed_receptor.pdb to Meeko
    mk_prepare_receptor.py -i fixed_receptor.pdb -o ${receptor_pdb.baseName} -p -a
    
    # 3. Prepare ligand
    mk_prepare_ligand.py -i ${ligand_sdf} -o ${ligand_sdf.baseName}.pdbqt
    """
}