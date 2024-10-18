//
// Create channel from folder
//

workflow INPUT_FROM_FOLDER {

    take:
    infolder

    main:
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
                num_fields_of_interest = "$params.sample_name_fields".toInteger()
                split_field_num = old_meta.split('_').size()
                num_fields_to_cut = split_field_num - num_fields_of_interest
                num_fields_to_cut = num_fields_to_cut == 0 ? 2 : num_fields_to_cut + 1
                sample_name = old_meta.split('_')[0..-num_fields_to_cut].join('_')
                meta.id = sample_name.replaceAll('-','_')
                meta.single_end = singleEnd
                lane_n_fastq = tuple((fastq.name =~ /L00\d/)[0], fastq)
                [meta, lane_n_fastq] }
        .groupTuple()
        .map{
            meta, lane_n_fastq ->
                meta = meta
                fastq = lane_n_fastq*.getAt(1).flatten()
                [meta, fastq] }

    emit:
    ch_fastq
}
