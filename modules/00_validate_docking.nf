process REDOCK_VALIDATION {
    label 'gpu_intensive'
    conda 'conda-forge::meeko conda-forge::pdb-tools'
    publishDir "${params.outdir}/00_validation", mode: 'copy'

    input:
    path wt_pdb        // estructura WT limpia (salida de CLEAN_STRUCTURE)
    path complex_pdb   // complejo cristalográfico original (para extraer la pose de referencia)
    path ligand_sdf    // AZ191.sdf (estructura del ligando)

    output:
    path "redocking_validation.txt", emit: report

    script:
    """
    # 1. Preparar receptor WT en PDBQT
    pdb_element ${wt_pdb} > fixed_wt.pdb
    mk_prepare_receptor.py -i fixed_wt.pdb -o wt_receptor -p -a

    # 2. Preparar ligando en PDBQT
    mk_prepare_ligand.py -i ${ligand_sdf} -o ligand.pdbqt

    # 3. Extraer el ligando cristalográfico del complejo (HETATM no-agua).
    #    Debe hacerse ANTES del docking: se usa como autobox_ligand para centrar
    #    la caja sobre el sitio de unión real y como referencia para el RMSD.
    grep -E "^HETATM" ${complex_pdb} | grep " ${params.ligand_resname} " > crystal_ligand.pdb
    echo "END" >> crystal_ligand.pdb

    # 4. Redocking de AZ191 sobre la WT.
    #    --autobox_ligand centra la caja en las coordenadas cristalográficas del ligando,
    #    evitando el problema de que AZ191.sdf libre esté centrado en el origen (≈0,0,0).
    gnina -r wt_receptor.pdbqt -l ligand.pdbqt \
        --autobox_ligand crystal_ligand.pdb \
        --autobox_add 4 \
        --cnn_scoring rescore \
        --seed ${params.seed} \
        --exhaustiveness ${params.exhaustiveness} \
        --num_modes ${params.num_modes} \
        --out redocked_wt.sdf

    # 5. Calcular RMSD y emitir informe (sale con error si RMSD > umbral)
    calc_rmsd.py \
        --ref-pdb  crystal_ligand.pdb \
        --docked   redocked_wt.sdf \
        --template ${ligand_sdf} \
        --threshold ${params.rmsd_threshold} \
        --output   redocking_validation.txt
    """
}
