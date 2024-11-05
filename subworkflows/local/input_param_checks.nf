// 
// Sanity checks on input parameters and creation of channels
// 

include { INPUT_FROM_FOLDER } from './input_from_folder.nf'
include { INPUT_FROM_SAMPLESHEET } from './input_from_samplesheet.nf'

workflow PARAMETER_CHECKS {

    take:
        ch_fasta
        ch_index
        ch_versions

    main:
        if (!params.bsgenome && (!params.forgeseed || !params.sourcedir)) {
            exit 1, 'Either the --bsgenome option or the following two options must be specified: --forgeseed, --sourcerdir.'
        } else if (params.bsgenome && (params.forgeseed || params.sourcedir)) {
            exit 1, 'The --bsgenome option and the following two options are mutually exclusive: --forgeseed, --sourcerdir.'
        }

        // if index is specified, it is used as input
        if (!params.fasta && !params.index) {
            exit 1, 'Reference FASTA file (--fasta) or genome index (--index) should be specified.'
        } else if (params.fasta && params.index) {
            exit 1, 'Only one of the two options, --fasta or --index, can be provided.'
        } else if (params.index) {
            // TODO: better id name?
            ch_index = Channel.of([
                        [ id:'index' ],
                        [file(params.index)]])
        } else {
            Channel
                .fromPath(params.fasta)
                .set { ch_fasta }
        }

        if (params.dist) {
            if (!params.dedup) {
                exit 1, 'The --dist option requires the --dedup option.'
            }
        }

        if (params.gtf) {
            // TODO: better id name?
            ch_gtf = Channel.of([
                    [ id:'gtf' ],
                    file(params.gtf, checkIfExists: true)])
        } else {
            exit 1, "The --gtf argument is mandatory."
        }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" & !params.fasta) {
            exit 1, 'The --splicesites option can only be used with the --fasta option.'
        }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" && params.bowtie2) {
            exit 1, 'The --splicesites option is mutually exclusive with the --bowtie2 option.'
        }

        if (params.samplesheet) {
            INPUT_FROM_SAMPLESHEET (
                params.samplesheet
            )
            ch_fastq = INPUT_FROM_SAMPLESHEET.out.ch_fastq
            ch_versions = ch_versions.mix(INPUT_FROM_SAMPLESHEET.out.versions)
        } else if (params.infolder) {
            ch_fastq = INPUT_FROM_FOLDER(
                params.infolder
            )
        } else {
            exit 1, 'Provide input by using the --samplesheet or the --infolder options.'
        }

    emit:
        ch_fasta
        ch_index
        ch_gtf
        ch_fastq
        ch_versions

}