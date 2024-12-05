// 
// Quality Control steps of CAGEr
// 

process CAGER_TAG_QC {
    label 'process_medium'
    stageInMode 'copy'

    input:
    path cager_obj
    path txdb
    path bsgenome_file
    val bsgenome_name

    output:
    path "annotated_cagexp.rds", emit: cager_rds
    path "corr_m.rds", emit: correlation_rds
    tuple path("*.pdf"), path("*plot.rds"), emit: plots
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cager_tag_qc.R  \
        -i ${cager_obj} \
        -a ${txdb} \
        -b \${bsgenome} \
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