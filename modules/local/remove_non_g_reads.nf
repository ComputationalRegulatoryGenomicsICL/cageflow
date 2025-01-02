process REMOVE_NON_G {
    tag "$meta.id"
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    tuple val(meta), path(fastqgz)

    output:
    path "*.fq.gz",      emit: reads
    path "versions.yml", emit: versions

    """
    fastqgzvar=${fastqgz}

    zcat ${fastqgz} | \\
    awk '{if (NR % 4 == 2 && substr(\$0, 1, 1) == "G") { \\
              print last; \\
              print \$0; \\
              getline; \\
              print \$0; \\
              getline; \\
              print \$0; \\
          } else { \\
              last=\$0; \\
          } \\
         }' | \\
    gzip > \${fastqgzvar%%.fq.gz}.with_g.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AWK: \$(awk -W version 2>&1 | head -1)
    END_VERSIONS
    """
}
