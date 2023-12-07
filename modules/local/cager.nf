process CAGER {
    label 'process_medium'
    stageInMode 'copy'

    // container 'docker://hub.docker.com/nikitinpavel/cager:0.5'
   
    input:
    path bsgnome
    val meta_bam

    output:
    path "*.RDS",        emit: rds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    echo ${meta_bam} | \\
        sed 's/, \\[/\\n/g' | \\
        tr -d '[],' | \\
        tr ' ' '\\t' | \\
        sed 's/id://' | \\
        sed 's/single_end://' \\
            > sample_list.tsv

    cager.R ${bsgnome} sample_list.tsv ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_BSgenome: \$(Rscript -e 'packageVersion("BSgenome")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}