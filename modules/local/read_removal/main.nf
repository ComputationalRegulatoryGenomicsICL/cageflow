//
// Removing reads which do not start with G
//

process READ_REMOVAL {
    label 'process_medium'
    stageInMode 'copy'

    input:
    tuple val(meta), path(input_fastq)

    output:
    tuple val(meta), path("filtered_output.fastq.gz")

    """
    remove_non_g_fastq.py -f ${input_fastq} -b 'G'
    """
}