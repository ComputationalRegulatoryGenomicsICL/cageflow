process CAGER_BAM {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    val bsgenome
    val meta_bam

    output:
    path "*.rds",        emit: rds
    path "versions.yml", emit: versions

    """
    echo ${meta_bam} | \\
        sed 's/, \\[/\\n/g' | \\
        tr -d '[],' | \\
        tr ' ' '\\t' | \\
        sed 's/id://' | \\
        sed 's/single_end://' \\
            > sample_list.tsv

    cager_bam.R ${bsgenome} sample_list.tsv ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_CAGEr: \$(Rscript -e 'packageVersion("CAGEr")' | awk '{print \$2}' | tr -d "‘’")
        R_BSgenome: \$(Rscript -e 'packageVersion("BSgenome")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}