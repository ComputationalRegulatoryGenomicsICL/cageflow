process BSGENOME {
    // tag "$meta.id"
    label 'process_low'
    stageInMode 'copy'

    // conda "YOUR-TOOLsrs/YOUR-TOOL-HERE' }"

    // input:
    // tuple val(meta), path(bam)

    output:
    path "*.tar.gz", emit: bsgenome
    // tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # http://hgdownload.cse.ucsc.edu/goldenPath/sacCer1/bigZips/sacCer1.fa.gz
    wget https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Scerevisiae.UCSC.sacCer1_1.4.0.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget 2>&1 | head -1 | cut -d" " -f1,2)
    END_VERSIONS
    """
}
