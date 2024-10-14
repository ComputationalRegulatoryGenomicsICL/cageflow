// 
// Sanity checks on input parameters and creation of channels
// 

include { INPUT_FROM_FOLDER } from './input_from_folder.nf'
include { INPUT_FROM_SAMPLESHEET } from './input_from_samplesheet.nf'

workflow PARAMETER_CHECKS {

    take:
        ch_fasta
        ch_index
        ch_gtf
        ch_txdb
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
            // if (params.gtf && params.txdb) {
            //     exit 1, 'When using STAR index as input, either --gtf or --txdb can be provided for CAGEr, but not both at the same time.'
            // }
            Channel
                .fromPath(params.index)
                .set { ch_index }
        } else {
            Channel
                .fromPath(params.fasta)
                .set { ch_fasta }
        }

        // either fasta or chromsizes file is required
        // Note: we do not check if the same fasta file is used to calculate the chromsizes and the index
        // if (!params.fasta && !params.chromsizes) {
        //     exit 1, 'Reference FASTA file (--fasta) or genome chromosome sizes (--chromsizes) should be specified.'
        // }

        if (params.dist) {
            if (!params.dedup) {
                exit 1, 'The --dist option requires the --dedup option.'
            }
        }

        if (params.gtf) {
            Channel
                .fromPath(params.gtf)
                .set { ch_gtf }
        } else {
            exit 1, "The --gtf argument is mandatory."
        }

        // if (params.gtf != "$projectDir/assets/NO_FILE_GTF" & !params.fasta) {
        //     exit 1, 'The --gtf option can only be used with the --fasta option.'
        // }

        // if (params.gtf != "$projectDir/assets/NO_FILE_GTF" & params.bowtie2) {
        //     exit 1, 'The --gtf option is mutually exclusive with the --bowtie2 option.'
        // }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" & !params.fasta) {
            exit 1, 'The --splicesites option can only be used with the --fasta option.'
        }

        if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" && params.bowtie2) {
            exit 1, 'The --splicesites option is mutually exclusive with the --bowtie2 option.'
        }

        if (params.chromsizes != "$projectDir/assets/NO_FILE_CHROMSIZES" & params.bowtie2) {
            exit 1, 'The --chromsizes option is mutually exclusive with the --bowtie2 option.'
        }

        // if (params.chromsizes == "$projectDir/assets/NO_FILE_CHROMSIZES" & !params.bowtie2) {
        //     exit 1, 'The use of the default mapper STAR requires the --chromsizes option.'
        // }

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
        ch_fastq
        ch_gtf
        ch_versions

}