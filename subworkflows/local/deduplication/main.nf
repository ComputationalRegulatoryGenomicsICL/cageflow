//
// Subworkflow for deduplication tasks
//

include { SAMTOOLS_SORT as SORT_FOR_FIXMATE} from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FIXMATE } from '../../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_SORT as SORT_AFTER_FIXMATE} from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as INDEX_AFTER_FIXMATE} from '../../modules/nf-core/samtools/index/main.nf'

include { SAMTOOLS_DEDUP } from '../../modules/local/samtools/dedup/main.nf'
include { SAMTOOLS_INDEX as INDEX_DEDUP} from '../../modules/nf-core/samtools/index/main.nf'

workflow DEDUPLICATION {
    take:
        ch_aligned
        ch_versions
        ch_for_cager

    main:
        SORT_FOR_FIXMATE (
            ch_aligned
        )
        ch_versions = ch_versions.mix(SORT_FOR_FIXMATE.out.versions)

        SAMTOOLS_FIXMATE (
            SORT_FOR_FIXMATE.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions)
        ch_bam_to_sort = SAMTOOLS_FIXMATE.out.bam

        SORT_AFTER_FIXMATE(ch_bam_to_sort)
        ch_versions = ch_versions.mix(SORT_AFTER_FIXMATE.out.versions)

        INDEX_AFTER_FIXMATE(SORT_AFTER_FIXMATE.out.bam)
        ch_versions = ch_versions.mix(
            INDEX_AFTER_FIXMATE.out.versions)

        SAMTOOLS_DEDUP (SORT_AFTER_FIXMATE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_DEDUP.out.versions)
        INDEX_DEDUP (SAMTOOLS_DEDUP.out.bam)
        ch_versions = ch_versions.mix(INDEX_DEDUP.out.versions)

        ch_bam_bai = SAMTOOLS_DEDUP.out.bam.join(INDEX_DEDUP.out.bai)

        if (params.bowtie2) {
            ch_for_cager = SAMTOOLS_DEDUP.out.bam
        }

    emit:
        ch_for_cager
        ch_bam_bai
        ch_versions
}