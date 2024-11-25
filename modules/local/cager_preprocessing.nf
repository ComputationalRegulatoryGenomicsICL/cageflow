// 
// Calling of tag clusters with CAGEr
// 

process CAGER_PREPROCESSING {
    label 'process_medium'
    stageInMode 'copy'

    input:
    path cager_obj

    output:
    path "*.rds",        emit: rds
    path "versions.yml", emit: versions

    """
    cager_preprocessing.R  \
        -i ${cager_obj} \
        -n ${params.norm_range_min} \
        -m ${params.norm_range_max} \
        -e ${params.norm_method} \
        -t ${params.total_tag_num} \
        -s ${params.sample_num_thr} \
        -r ${params.ctss_thr} \
        -p ${projectDir} \
        -c ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_packages: \$(Rscript -e 'sessionInfo(package = NULL)')
    END_VERSIONS
    """
}