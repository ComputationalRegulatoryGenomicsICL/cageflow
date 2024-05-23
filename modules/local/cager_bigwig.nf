process CAGER_BIGWIG {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    val bsgenome
    val bigwig

    output:
    path "*.rds",        emit: rds
    path "versions.yml", emit: versions

    """
    cager_bigwig.R ${bsgenome} "${bigwig}" ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_BSgenome: \$(Rscript -e 'packageVersion("BSgenome")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}