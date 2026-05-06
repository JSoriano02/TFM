process FILTER_MUTATIONS {
    // Instruct Nextflow to use Conda to install pandas automatically
    conda 'conda-forge::pandas'

    input:
    path gdc_tsv

    output:
    path "filtered_mutations.csv", emit: mutations_csv

    script:
    """
    # Nextflow automatically finds filter_mutations.py in the bin/ folder
    filter_mutations.py -i ${gdc_tsv} -o filtered_mutations.csv -n 3
    """
}