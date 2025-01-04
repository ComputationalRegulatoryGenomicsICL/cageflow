// 
// Processing before mapping
// 

include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../../modules/nf-core/fastqc/main.nf'
include { TRIMGALORE } from '../../modules/nf-core/trimgalore/main.nf'
include { REMOVE_NON_G } from '../../modules/local/remove_non_g_reads.nf'
include { CUTADAPT } from '../../modules/nf-core/cutadapt/main.nf'

workflow PREPROCESSING {

    take:
        ch_fastq
        ch_versions
        ch_multiqc_files

    main:
        CAT_FASTQ (
            ch_fastq
        ).reads.set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

        FASTQC (
            ch_cat_fastq
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions)

        TRIMGALORE (
            ch_cat_fastq
        )
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

        if (params.removenong) {
            REMOVE_NON_G (
                TRIMGALORE.out.reads
            )
            ch_versions = ch_versions.mix(REMOVE_NON_G.out.versions)
        }

        if (!params.nogtrim) {
            if (!params.removenong) {
                CUTADAPT (
                    TRIMGALORE.out.reads
                )
            } else {
                CUTADAPT (
                    REMOVE_NON_G.out.reads
                )
            }
            ch_versions = ch_versions.mix(CUTADAPT.out.versions)
        }

        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))

        if (!params.nogtrim) {
            ch_reads_to_align = CUTADAPT.out.reads
        } else if (params.removenong) {
            ch_reads_to_align = REMOVE_NON_G.out.reads
        } else {
            ch_reads_to_align = TRIMGALORE.out.reads
        }

        //ch_reads_to_align = !params.nogtrim ? CUTADAPT.out.reads : TRIMGALORE.out.reads

    emit:
        ch_reads_to_align
        ch_versions
        ch_multiqc_files

}