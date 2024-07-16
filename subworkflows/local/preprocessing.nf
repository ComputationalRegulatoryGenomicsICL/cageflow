// 
// Processing before mapping
// 
include { FORGE_BSGENOME } from '../modules/local/forge_bsgenome.nf'

include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../modules/nf-core/fastqc/main.nf'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
include { CUTADAPT } from '../modules/nf-core/cutadapt/main.nf'