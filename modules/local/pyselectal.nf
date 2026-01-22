process PYSELECTAL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pysam=0.23.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_bam = "${prefix}.pyselectal.bam"

    def args = []
    if (params.pyselectal?.min_softclip != null) args << "-n ${params.pyselectal.min_softclip}"
    if (params.pyselectal?.max_softclip != null) args << "-m ${params.pyselectal.max_softclip}"
    if (params.pyselectal?.prefix)               args << "-p '${params.pyselectal.prefix}'"
    if (params.pyselectal?.k)                    args << "-k ${params.pyselectal.k}"
    if (params.pyselectal?.paired)               args << "--paired"
    if (params.pyselectal?.sort)                 args << "--sort"
    def threads = params.pyselectal?.threads ?: task.cpus
    if (threads) args << "-t ${threads}"

    def extra = task.ext.args ?: ''
    def cli = (args + [extra]).findAll{ it && it.toString().trim() }.join(' ')

    
    """
    pyselectal.py \\
        ${cli} \\
        ${bam} \\
        ${out_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: "$(python -c 'import pysam; print(getattr(pysam, "__version__", "1.0"))')"
    END_VERSIONS
    """
}
