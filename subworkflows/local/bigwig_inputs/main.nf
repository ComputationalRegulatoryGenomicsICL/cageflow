//
// Create channel from mapping file content
//

workflow BIGWIG_INPUTS {
    take:
    sample_file

    main:
    input_files = sample_file
        .splitCsv( header:true , sep:',')
        .map { create_sample_channel(it) }

    emit:
    input_files
}

def create_sample_channel(LinkedHashMap row) {
    bigwigs = row.path
    bigwig_1 = file(bigwigs.split(' ')[0].minus('['))
    bigwig_2 = file(bigwigs.split(' ')[1].minus(']'))

    return [bigwig_1, bigwig_2]
}