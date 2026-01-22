include { PYSELECTAL } from '../../modules/local/pyselectal.nf'
// include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main.nf

workflow FILTER_BAM_OPTIONAL {

    take:
        ch_aligned
        ch_versions

    main:
        if (params.pyselectal?.enabled) {
            PYSELECTAL(ch_aligned)
            ch_versions = ch_versions.mix(PYSELECTAL.out.versions)
            ch_aligned_out = PYSELECTAL.out.bam
        } else {
            ch_aligned_out = ch_aligned
        }

    emit:
        ch_aligned = ch_aligned_out
        ch_versions = ch_versions
}