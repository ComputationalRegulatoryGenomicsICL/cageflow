// 
// Quality Control steps of CAGEr
// 

process CAGER_TAG_QC {
    label 'process_verylong'
    stageInMode 'copy'

    input:
    path cager_obj
    path txdb
    path bsgenome_file
    val bsgenome_name

    output:
    path "intermediate_cagerobj/annotated_cagexp.rds", emit: cager_rds
    path "plots/corr_m.rds", emit: correlation_rds
    tuple path("plots/*.pdf"), path("plots/*plot.rds"), emit: plots
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cager_tag_qc.R  \
        --cageexp_object ${cager_obj} \
        --annotation ${txdb} \
        --bsgenome \${bsgenome} \
        --corrplot_tagCountThreshold ${params.corrplot_tagCountThreshold} \
        --heatmap_cex ${params.heatmap_cex} \
        --project_dir ${projectDir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}