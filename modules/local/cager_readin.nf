// Read in to CAGEr in bigwig or bam format

process CAGER_READIN {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    path bsgenome_file
    val bsgenome_name
    val data_in
    val data_type

    output:
    path "*.rds",        emit: rds
    path "versions.yml", emit: versions

    """
    if [ -z ${bsgenome_name} ]
    then
        bsgenome=${bsgenome_file}
    else
        bsgenome=${bsgenome_name}
    fi

    echo ${data_in} | \\
        sed 's/, \\[/\\n/g' | \\
        tr -d '[],' | \\
        tr ' ' '\\t' | \\
        sed 's/id://' | \\
        sed 's/single_end://' \\
            > sample_list.tsv

    cager_readin.R \
        -t "${data_type}" \
        -b \${bsgenome} \
        -s sample_list.tsv \
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