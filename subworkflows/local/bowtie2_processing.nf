// 
// Subworkflow to define flow from Bowtie2
// 

include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main.nf' 
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main.nf'
include { SAMTOOLS_VIEW_MAPQ } from '../modules/nf-core/samtools/view_mapq/main.nf'
include { SAMTOOLS_FIXMATE } from '../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_DEDUP } from '../modules/local/samtools_dedup.nf'
include { CAGER_BAM } from '../modules/local/cager_bam.nf'