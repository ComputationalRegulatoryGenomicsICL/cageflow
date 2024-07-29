// 
// Quality Control steps of CAGEr
// 

process CAGER_TAG_QC {
    label 'process_medium'
    stageInMode 'copy'

    input:
    path cager_obj

    output:
    path "versions.yml", emit: versions

    """
    cager_tag_qc.R  \
        -i ${cager_obj} \
        -a ${params.annotation} \
        -p ${projectDir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}