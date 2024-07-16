// 
// Subworkflow to define flow from STAR
// 

include { STAR_ALIGN } from '../modules/nf-core/star/align/main.nf' 
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main.nf'
include { UCSC_WIGTOBIGWIG } from '../modules/nf-core/ucsc/wigtobigwig/main.nf' 
include { CAGER_BIGWIG } from '../modules/local/cager_bigwig.nf'