process FOLDX_MUTAGENESIS {
    conda 'conda-forge::python=3.10 conda-forge::pandas'
    publishDir "${params.outdir}/03_mutants", mode: 'copy'

    input:
    path wt_pdb
    path mutations_csv

    output:
    path "*_mutant.pdb",          emit: mutant_pdbs
    path "*_ddg_replicas.csv",    emit: ddg_replicas
    path "foldx_ddg_summary.csv", emit: ddg_summary

    script:
    """
    # 1. Create a physical copy while simultaneously removing any Windows
    # carriage returns (CRLF to LF) that crash FoldX
    tr -d '\\r' < ${wt_pdb} > ./physical_WT.pdb

    # 2. Force read-write permissions (FoldX C++ core bug workaround)
    chmod 644 ./physical_WT.pdb

    # 3. Check if the file was processed successfully and is not empty
    if [ ! -s physical_WT.pdb ]; then
        echo "Critical Error: physical_WT.pdb is empty. Cleaning step may have failed."
        exit 1
    fi

    # 4. Verify and copy rotabase.txt
    ROTABASE_PATH="/home/wormi/software/rotabase.txt"
    if [ ! -f "\$ROTABASE_PATH" ]; then
        echo "Critical Error: rotabase.txt not found in \${ROTABASE_PATH}"
        exit 1
    fi
    cp "\$ROTABASE_PATH" ./rotabase.txt

    # 5. Execute the python script passing the fixed physical copy and number of replicas
    run_foldx.py -i ${mutations_csv} -p physical_WT.pdb --runs ${params.foldx_runs}
    """
}