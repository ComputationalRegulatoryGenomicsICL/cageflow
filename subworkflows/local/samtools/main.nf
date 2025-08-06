//
// Processing of files after mapping
//

include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'

include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main.nf'

workflow SAMTOOLS {
    take:
        ch_aligned
        ch_versions
        ch_for_cager

    main:
        SAMTOOLS_SORT(ch_aligned)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
        SAMTOOLS_INDEX (SAMTOOLS_SORT.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        if (params.bowtie2) {
            ch_for_cager = SAMTOOLS_SORT.out.bam
        }

    emit:
        ch_for_cager
        ch_bam_bai
        ch_versions

}
