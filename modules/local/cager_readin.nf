// Read in to CAGEr in bigwig or bam format

process CAGER_READIN {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    path bsgenome_file
    val bsgenome_name
    val sample_table
    val data_type
    path ch_collected

    output:
    path "intermediate_cagerobj/*.rds",        emit: rds
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    cager_readin.R \
        -t "${data_type}" \
        -b \${bsgenome} \
        -s "${sample_table}" \
        -p ${projectDir} \
        -c ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_BSgenome: \$(Rscript -e 'packageVersion("BSgenome")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}