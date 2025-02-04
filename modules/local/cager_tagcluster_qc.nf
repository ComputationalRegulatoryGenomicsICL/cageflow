// 
// Quality Control steps of CAGEr after clustering of CTSS
// 

process CAGER_TAGCLUSTER_QC {
    label 'process_medium'
    stageInMode 'copy'

    input:
    path cager_obj
    path txdb
    path bsgenome_file
    val bsgenome_name

    output:
    path "tables/*.csv", emit: counts_csv
    tuple path("plots/*.pdf"), path("plots/*plot.rds"), emit: plots
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cager_tagcluster_qc.R  \
        -i ${cager_obj} \
        -a ${txdb} \
        -b \${bsgenome} \
        -o ${params.iq_low} \
        -g ${params.iq_high} \
        -u ${params.tssregion_up} \
        -d ${params.tssregion_down} \
        -l ${params.tsslogo_upstream} \
        -p ${projectDir} \
        -k ${params.pca_rank}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}