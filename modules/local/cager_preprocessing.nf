// 
// Calling of tag clusters with CAGEr
// 

process CAGER_PREPROCESSING {
    label 'process_verylong'
    stageInMode 'copy'

    input:
    path cager_obj
    path bsgenome_file
    val bsgenome_name
    path txdb

    output:
    path "intermediate_cagerobj/normalized_clustered_cagexp.rds",        emit: rds
    tuple path("plots/*.pdf"), path("plots/*.txt"), path("plots/*plot.rds"), emit: results
    tuple path("tracks/*.bw"), path("tracks/*.bed"), emit: tracks
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cager_preprocessing.R  \
        -i ${cager_obj} \
        -n ${params.norm_range_min} \
        -m ${params.norm_range_max} \
        -e ${params.norm_method} \
        -t ${params.T_norm} \
        -s ${params.sample_num_thr} \
        -r ${params.ctss_thr} \
        -u ${params.consensus_ctss_thr} \
        -d ${params.consensus_ctss_dist} \
        -a ${txdb} \
        -p ${projectDir} \
        -b \${bsgenome} \
        -c ${task.cpus}
    
    rm Rplots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}