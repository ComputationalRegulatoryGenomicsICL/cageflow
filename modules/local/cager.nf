process CAGER {
    tag "$meta.id"
    label 'process_medium'
    stageInMode 'copy'

    container 'nikitinpavel/cager:0.1'

    input:
    tuple val(meta), path(bam)
    path rscript

    output:
    tuple val(meta), path("*.RDS"), emit: rds
    // tuple val(meta), path("tsv")   , emit: tsv_dir
    // tuple val(meta), path("pdf")   , emit: pdf_dir
    // tuple val(meta), path("*.html"), emit: knitted_html
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cager.R ${bam}

    # render_rmd.R $rmd bowtie2_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}
