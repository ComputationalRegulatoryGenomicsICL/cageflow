// 
// Calling of enhancers with CAGEfightr
// 

process CAGEFIGHTR_ENHANCERS {
    label 'process_verylong'
    stageInMode 'copy'

    input:
    path cager_obj
    path txdb
    path correlation_rds
    path bsgenome_file
    val bsgenome_name

    output:
    tuple path("intermediate_cagerobj/supported_enhancers.rds"), path("intermediate_cagerobj/nonTSS_enhancers.rds"), emit: rds
    tuple path("plots/*.pdf"), path("plots/*plot.rds"), emit: plots
    tuple path("tables/*.tsv"), path("tracks/*.bed"), emit: enhancer_table
    path "versions.yml", emit: versions
    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cagefightr_enhancer_calling.R  \
        --cageexp_object ${cager_obj} \
        --annotation ${txdb} \
        --cfBalanceThreshold ${params.cfBalanceThreshold} \
        --unexpressed ${params.unexpressed} \
        --minSamples ${params.minSamples} \
        --remove_gg_initiator ${params.remove_gg_initiator} \
        --keep_only_yr_yc ${params.keep_only_yr_yc} \
        --tssregion_up ${params.tssregion_up} \
        --tssregion_down ${params.tssregion_down} \
        --bsgenome \${bsgenome} \
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