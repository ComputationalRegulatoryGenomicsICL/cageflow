// 
// Subworkflow to define flow from STAR
// 

include { STAR_ALIGN } from '../../modules/nf-core/star/align/main.nf' 
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main.nf'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main.nf' 

workflow STAR_PROCESSING {

    take:
        ch_reads_to_align
        ch_fasta
        ch_index
        ch_chrom_sizes
        ch_multiqc_files
        ch_versions

    main:
        if (!params.index) {
            gtf_file = file(params.gtf, checkIfExists: true)
            splice_sites_file = file(params.splicesites, checkIfExists: true)
            STAR_GENOMEGENERATE (
                ch_fasta,
                gtf_file,
                splice_sites_file
            )
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            
            ch_index = STAR_GENOMEGENERATE.out.index
        }

        STAR_ALIGN (
            ch_reads_to_align,
            ch_index
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        ch_aligned = STAR_ALIGN.out.bam_sorted

        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

        UCSC_WIGTOBIGWIG (
            STAR_ALIGN.out.wigtobigwig,
            ch_chrom_sizes
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions)

        bigwig_ch_for_cager = UCSC_WIGTOBIGWIG.out.bw
            .map { it[1] }
            .collect()

    emit:
        bigwig_ch_for_cager
        ch_aligned
        ch_multiqc_files
        ch_versions


}