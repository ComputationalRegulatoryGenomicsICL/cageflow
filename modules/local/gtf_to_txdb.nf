process GTF_TO_TXDB {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    path gtf_file

    output:
    path "*.sqlite",     emit: txdb
    path "versions.yml", emit: versions

    """
    gtf_to_txdb.R ${gtf_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_BSgenome: \$(Rscript -e 'packageVersion("txdbmaker")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}