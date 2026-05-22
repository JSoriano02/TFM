process CLEAN_STRUCTURE {
    conda 'conda-forge::python=3.10'
    publishDir "${params.outdir}/02_cleaned", mode: 'copy'

    input:
    path raw_pdb

    output:
    path "${raw_pdb.baseName}_clean.pdb", emit: cleaned_pdb

    script:
    """
    # Use the custom robust Python cleaner for FoldX 4 compatibility
    clean_pdb_foldx.py -i ${raw_pdb} -o ${raw_pdb.baseName}_clean.pdb
    """
}