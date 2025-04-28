// 
// Calling of enhancers with CAGEfightr
// 

process CAGEFIGHTR_ENHANCERS {
    label 'process_verylong'
    stageInMode 'copy'

    input:
    path cager_obj
    path txdb

    output:
    path "intermediate_cagerobj/enhancers.rds",        emit: rds
    tuple path("plots/*.pdf"), path("plots/*plot.rds"), emit: plots
    path "tables/*.tsv", emit: enhancer_table
    path "versions.yml", emit: versions
    """

    cagefightr_enhancer_calling.R  \
        --cageexp_object ${cager_obj} \
        --annotation ${txdb} \
        --cfBalanceThreshold ${params.cfBalanceThreshold}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}