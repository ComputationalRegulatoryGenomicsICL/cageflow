//
// Sanity checks on input parameters and creation of channels
//

include { samplesheetToList         } from 'plugin/nf-schema'
include { INPUT_FROM_FOLDER         } from '../input_from_folder/main.nf'

workflow PARAMETER_CHECKS {

    take:
        ch_fasta
        ch_index
        ch_versions

    main:

        if (params.samplesheet) {
            //
            // Create channel from samplesheet file provided through params.samplesheet
            //

            Channel
                .fromList(samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json"))
                .map {
                    meta, fastq_1, fastq_2 ->
                        if (!fastq_2) {
                            return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                        } else {
                            return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                        }
                }
                .groupTuple()
                .map { samplesheet ->
                    validateInputSamplesheet(samplesheet)
                }
                .map {
                    meta, fastq ->
                        meta.id = meta.id.split('_')[0..-2].join('_')
                        [ meta, fastq ] }
                    .groupTuple(by: [0])
                    .map{ meta, fastq -> [ meta, fastq.flatten() ] 
                }
                .set { ch_fastq }

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
        if (!params.genome && !params.index) {
            exit 1, 'Reference genome FASTA file (--fasta) or genome index (--index) should be specified.'
        } else if (params.index) {
            ch_pre_idx = Channel.fromPath(params.index, checkIfExists: true)
            ch_index = sample_meta.combine(ch_pre_idx)
            if (params.genome) {
                ch_pre_fa = Channel.fromPath(params.genome, checkIfExists: true)
                ch_fasta = ch_genome_name.combine(ch_pre_fa)
            } else {
                ch_fasta = Channel.empty()
            }
        } else {
            ch_pre_fa = Channel.fromPath(params.genome, checkIfExists: true)
            ch_fasta = ch_genome_name.combine(ch_pre_fa)
            ch_index = Channel.empty()
        }

        if (params.dist) {
            if (!params.dedup) {
                exit 1, 'The --dist option requires the --dedup option.'
            }
        }

        if (!params.bsgenome && (!params.forgeseed || !params.sourcedir)) {
            exit 1, 'Either the --bsgenome option or the following two options must be specified: --forgeseed, --sourcerdir.'
        } else if (params.bsgenome && (params.forgeseed || params.sourcedir)) {
            exit 1, 'The --bsgenome option and the following two options are mutually exclusive: --forgeseed, --sourcerdir.'
        }


    emit:
        ch_fasta
        ch_index
        ch_fastq
        ch_versions

}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
