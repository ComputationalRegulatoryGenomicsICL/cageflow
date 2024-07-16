// 
// Steps for creating summary statistics for quality check
// 


include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as SORT_FOR_FIXMATE} from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP} from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_STATS } from '../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_IDXSTATS } from '../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/samtools/flagstat/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'