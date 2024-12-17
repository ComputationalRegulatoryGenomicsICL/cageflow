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

        sample_meta = ch_fastq.map{ meta, fastq ->
            meta = meta
            [meta]}

        ch_genome_name = Channel.of(params.genome_name)

        // // if index is specified, it is used as input
        if (!params.fasta && !params.index) {
            exit 1, 'Reference FASTA file (--fasta) or genome index (--index) should be specified.'
        } else if (params.index) {
            ch_pre_idx = Channel.fromPath(params.index, checkIfExists: true)
            ch_index = sample_meta.combine(ch_pre_idx)
            // NOTE: this will mean bowtie2 does not run
            // TODO: mock some input maybe, fasta is only needed for a specific type of output which we don't use
            if (params.fasta) {
                ch_pre_fa = Channel.fromPath(params.fasta, checkIfExists: true)
                ch_fasta = ch_genome_name.combine(ch_pre_fa)
            } else {
                ch_fasta = Channel.empty()
            }
        } else {
            //ch_fasta = Channel.fromPath(params.fasta)
            ch_pre_fa = Channel.fromPath(params.fasta, checkIfExists: true)
            ch_fasta = ch_genome_name.combine(ch_pre_fa)
            ch_index = Channel.empty()
        }

        // All hell broke loose with removing this option, so that now it is placed back
        // ch_pre_fa = Channel.fromPath(params.fasta)
        // ch_fasta = sample_meta.combine(ch_pre_fa)
        // ch_pre_idx = Channel.fromPath(params.index)
        // ch_index = sample_meta.combine(ch_pre_idx)

        if (params.gtf) {
            gtf_path = file(params.gtf)
            if (gtf_path.extension == 'gz') {
                ch_gtf = 
            } else {
                ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
            }
        } else {
            exit 1, "The --gtf argument is mandatory."
        }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" & !params.fasta) {
            exit 1, 'The --splicesites option can only be used with the --fasta option.'
        }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" && params.bowtie2) {
            exit 1, 'The --splicesites option is mutually exclusive with the --bowtie2 option.'
        }

        if (!params.bsgenome && (!params.forgeseed || !params.sourcedir)) {
            exit 1, 'Either the --bsgenome option or the following two options must be specified: --forgeseed, --sourcerdir.'
        } else if (params.bsgenome && (params.forgeseed || params.sourcedir)) {
            exit 1, 'The --bsgenome option and the following two options are mutually exclusive: --forgeseed, --sourcerdir.'
        }


    emit:
        ch_fasta
        ch_index
        ch_gtf
        ch_fastq
        ch_versions

}