// 
// Subworkflow to get the BSgenome via forging or loading
// 

include { FORGE_BSGENOME } from '../../modules/local/forge_bsgenome.nf'
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main.nf'

workflow PREPARE_METADATA {

    main:
        // prepare or fetch bsgenome 
        if (params.forgeseed) {
            forge_seed = file(params.forgeseed, checkIfExists: true)
            seqs_srcdir = file(params.sourcedir, checkIfExists: true)
            FORGE_BSGENOME (
                forge_seed,
                seqs_srcdir
            )
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

        if (!params.chromsizes){
            // prepare chromosome sizes
            chrom_name_ch = Channel.of(params.chrom_name)
            chrom_ch = Channel.fromPath(params.fasta)
            chromsize_ch = chrom_name_ch.combine(chrom_ch)

            CUSTOM_GETCHROMSIZES( chromsize_ch )

            ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
                .map{
                    name, size_file ->
                    size_file
                }
            versions = CUSTOM_GETCHROMSIZES.out.versions
        }
        

    emit:
        ch_bsgenome_file
        ch_bsgenome_name
        ch_chrom_sizes
        versions
}