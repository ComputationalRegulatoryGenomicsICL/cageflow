// 
// Write sample list to file
// 

process WRITE_SAMPLE_LIST {

    input:
    tuple val(meta), path(bw_or_bam)

    output:
    path("sample_list.tsv")

    script:
    if ( bw_or_bam[1] != null )
        """
        line="${meta.id},${meta.single_end},[${PWD}/${bw_or_bam[0]} ${PWD}/${bw_or_bam[1]}]" 
        echo \$line > sample_list.tsv
        """
    else
        """
        line="${meta.id},${meta.single_end},[${PWD}/${bw_or_bam[0]}]"
        echo \$line > sample_list.tsv
        """

}