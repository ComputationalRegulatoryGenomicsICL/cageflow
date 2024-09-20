// 
// Subworkflow to get the BSgenome via forging or loading
// 

include { FORGE_BSGENOME } from '../../modules/local/forge_bsgenome.nf'
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main.nf'

workflow PREPARE_METADATA {

    take:
        ch_versions

    main:

        // prepare or fetch bsgenome 
        if (params.forgeseed) {
            forge_seed = file(params.forgeseed, checkIfExists: true)
            seqs_srcdir = file(params.sourcedir, checkIfExists: true)
            FORGE_BSGENOME (
                forge_seed,
                seqs_srcdir
            )
            ch_versions = ch_versions.mix(FORGE_BSGENOME.out.versions)
        }

        if (params.bsgenome) {
            if (params.bsgenome.endsWith('.tar.gz')) {
                ch_bsgenome_file = file(
                    params.bsgenome,
                    checkIfExists: true)
                ch_bsgenome_name = ''
            } else {
                ch_bsgenome_file = file(
                    "$projectDir/assets/NO_FILE_BSGENOME")
                ch_bsgenome_name = params.bsgenome
            }
        } else {
            ch_bsgenome_file = FORGE_BSGENOME.out.bsgenome
            ch_bsgenome_name = ''
        }

        // prepare chromosome sizes
        if (!params.chromsizes){
            if (params.fasta) {
                chrom_ch = Channel.fromPath(params.fasta)

                CUSTOM_GETCHROMSIZES( chrom_ch )
                ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes

                ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
            } else { // a genome index was provided instead
                ch_chrom_sizes = Channel.fromPath(params.index + '/chrNameLength.txt')
            }
        } else {
            ch_chrom_sizes = Channel.fromPath(params.chromsizes)
        }

    emit:
        ch_bsgenome_file
        ch_bsgenome_name
        ch_chrom_sizes
        ch_versions
}