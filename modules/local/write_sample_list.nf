// 
// Write sample list to file
// 

process WRITE_SAMPLE_LIST {

    input:
    tuple val(meta), path(bw_or_bam)

    output:
    path("sample_list.tsv")

    shell:
    '''
    line="!{meta.id},!{meta.single_end},[!{bw_or_bam[0]} !{bw_or_bam[1]}]"
    echo $line >> sample_list.tsv
    '''

}