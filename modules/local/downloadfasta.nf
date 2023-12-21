process DOWNLOAD_FASTA {
    label 'process_single'
    stageInMode 'copy'

    input:
    val ucscid

    output:
    path "*.fa.gz",      emit: fasta
    path "versions.yml", emit: versions
    
    """
    wget "http://hgdownload.cse.ucsc.edu/goldenPath/${ucscid}/bigZips/${ucscid}.fa.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget 2>&1 | head -1 | cut -d" " -f1,2)
    END_VERSIONS
    """
}
