nextflow.enable.dsl=2

// Import the filtering module
include { FILTER_MUTATIONS } from './modules/01_filter_data.nf'

// Define input data parameters using the raw_dir from nextflow.config
params.mutations_tsv = "${params.raw_dir}/DYRK1B_mutations.tsv"

workflow {
    // Create a channel that emits the downloaded TSV file
    gdc_data_ch = Channel.fromPath(params.mutations_tsv)

    // Run the filtering process
    filtered_data = FILTER_MUTATIONS(gdc_data_ch)
    
    // Print a confirmation message to the terminal when finished
    filtered_data.mutations_csv.view { "The filtered file is ready at: $it" }
}