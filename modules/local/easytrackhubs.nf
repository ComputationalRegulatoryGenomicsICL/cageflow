process EASY_TRACK_HUBS {
    label 'process_single'
    stageInMode 'copy'

    input:
    path(normalized_bw)
    val(ref_genome)

    output:
    path "question", emit: trackhubs
    
    """
    make_trackhubs.R --bigwigs ${normalized_bw} --ref_genome ${ref_genome} 
    """
}