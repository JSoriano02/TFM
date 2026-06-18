// main.nf
nextflow.enable.dsl=2

include { REDOCK_VALIDATION } from './modules/00_validate_docking.nf'
include { FILTER_MUTATIONS }  from './modules/01_filter_data.nf'
include { CLEAN_STRUCTURE }   from './modules/02_clean_structure.nf'
include { FOLDX_MUTAGENESIS } from './modules/03_mutagenesis.nf'
include { PREPARE_MEEKO }     from './modules/04_prepare_meeko.nf'
include { DOCKING_GNINA }     from './modules/05_docking_gnina.nf'
include { EXTRACT_AND_REPORT }from './modules/06_analysis.nf'

params.mutations_tsv = "${params.raw_dir}/DYRK1B_mutations.tsv"

workflow {
    // 1. Filter Top 3 Mutations
    filtered_data = FILTER_MUTATIONS(file(params.mutations_tsv))
    
    // 2. Clean WT Protein
    cleaned_wt = CLEAN_STRUCTURE(file(params.wt_pdb))

    // 0 (control). Validar el protocolo de docking redocking AZ191 sobre la WT cristalográfica
    // Si RMSD > params.rmsd_threshold el proceso falla y el pipeline se detiene aquí.
    validation_out = REDOCK_VALIDATION(cleaned_wt, file(params.wt_pdb), file(params.ligand))

    // Enganchar la mutagénesis a la validación: FOLDX_MUTAGENESIS no arranca hasta que
    // REDOCK_VALIDATION haya completado con éxito.
    gated_wt  = cleaned_wt.combine(validation_out.report).map { wt, _r -> wt }

    // 3. Mutagenesis with FoldX (Outputs mutant PDBs + réplicas ΔΔG)
    foldx_out     = FOLDX_MUTAGENESIS(gated_wt, filtered_data)

    // Combine WT and Mutants into a single channel for preparation
    all_receptors = cleaned_wt.concat(foldx_out.mutant_pdbs.flatten())
    
    // 4. Prepare Ligand and Receptors with Meeko
    prepared = PREPARE_MEEKO(all_receptors, file(params.ligand))
    
    // 5. Docking with Gnina
    docking_results = DOCKING_GNINA(prepared.receptor_pdbqt, prepared.ligand_pdbqt)
    
    // 6. Informe completo: afinidad GNINA, estabilidad FoldX, interacciones ProLIF,
    //    clasificación y estadística con barras de error
    EXTRACT_AND_REPORT(docking_results.collect(), cleaned_wt,
                       foldx_out.mutant_pdbs.collect(), foldx_out.ddg_summary,
                       foldx_out.ddg_replicas.collect())
}