//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main.nf'

workflow INPUT_FROM_SAMPLESHEET {
    take:
    samplesheet // /path/to/samplesheet.csv

    main:
    input_handler = file(samplesheet, checkIfExists: true)
    reads = SAMPLESHEET_CHECK ( input_handler )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }

    ch_fastq = reads
        .map {
                meta, fastq ->
                    meta.id = meta.id.split('_')[0..-2].join('_')
                    [ meta, fastq ] }
            .groupTuple(by: [0])
            .map{ meta, fastq -> [ meta, fastq.flatten() ] }

    emit:
    ch_fastq                                  // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Format: [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []

    if ((file(row.fastq_1) == []) || (!file(row.fastq_1).exists())) {
        throw new Exception("Please check input samplesheet: Read 1 FastQ file does not exist!\nRead 1 FastQ file: ${row.fastq_1}")
    }

    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if ((file(row.fastq_2) == []) || (!file(row.fastq_2).exists())) {
            throw new Exception("Please check input samplesheet: Read 2 FastQ file does not exist!\nRead 2 FastQ file: ${row.fastq_2}")
        }

        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }

    return fastq_meta
}
