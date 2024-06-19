process FORGE_BSGENOME {
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    path forge_seed
    path seqs_srcdir

    output:
    path "BSgenome.*",   emit: bsgenome
    path "versions.yml", emit: versions

    """
    forge_bsgenome.R ${forge_seed} ${seqs_srcdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_BSgenome: \$(Rscript -e 'packageVersion("BSgenome")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
    """
}