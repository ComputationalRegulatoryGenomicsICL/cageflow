//
// Create channel from folder
//

workflow INPUT_FROM_FOLDER {

    take:
    infolder

    main:
    // check here if singleEndness is correct
    any_R2_file = file("$infolder/**_R2*fastq.gz")
    singleEnd = true
    if (any_R2_file.size() > 0){
        singleEnd = false
    }

    ch_fastq = channel
        .fromFilePairs(
            "$infolder/**_R{1,2}*fastq.gz",
            size: singleEnd ? 1 : 2)
        .map{
            old_meta, fastq -> 
                def meta = [:]
                meta.id = old_meta.split('_')[0..-3].join('_')
                meta.single_end = singleEnd
                fastq = tuple((fastq.name =~ /L00\d/)[0], fastq)
                [meta, fastq ] }
        .groupTuple()
        .map{
            meta, lane_n_fastq ->
                meta = meta
                fastq = lane_n_fastq*.getAt(1).flatten()
                [meta, fastq] }

    emit:
    ch_fastq
}
