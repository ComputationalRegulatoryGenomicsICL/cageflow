process REMOVE_NON_G {
    tag "$meta.id"
    label 'process_medium'
    stageInMode 'copy'
   
    input:
    tuple val(meta), path(fastqgz)

    output:
    tuple val(meta), path ("*.fq.gz"), emit: reads
    path "versions.yml",               emit: versions

    script:
    if (meta.single_end) {
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
    } else {
        """
        fastqgzvar0=${fastqgz[0]}
        fastqgzvar1=${fastqgz[1]}

        paste <(zcat ${fastqgz[0]}) <(zcat ${fastqgz[1]}) | \\
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
             }' > \\
        tmp.with_g.tsv

        <tmp.with_g.tsv cut -d\$'\\t' -f1 > \\
        \${fastqgzvar0%%.fq.gz}.with_g.fq

        <tmp.with_g.tsv cut -d\$'\\t' -f2 > \\
        \${fastqgzvar1%%.fq.gz}.with_g.fq

        gzip \${fastqgzvar0%%.fq.gz}.with_g.fq \${fastqgzvar1%%.fq.gz}.with_g.fq

        rm tmp.with_g.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            AWK: \$(awk -W version 2>&1 | head -1)
        END_VERSIONS
        """
    }
}
