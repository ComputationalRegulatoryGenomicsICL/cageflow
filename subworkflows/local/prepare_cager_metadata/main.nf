//
// Subworkflow to get the BSgenome via forging or loading
//

include { GTF2TXDB } from '../../../modules/local/gtf2txdb/main.nf'
include { FORGE_BSGENOME } from '../../../modules/local/forge_bsgenome/main.nf'

workflow PREPARE_CAGER_METADATA {

    take:
        ch_gtf
        ch_versions

    main:

        println("Prepare CAGEr metadata")

        // prepare or fetch BSgenome
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

        ch_txdb = GTF2TXDB(ch_gtf)
        ch_txdb_file = GTF2TXDB.out.txdb
        ch_versions = ch_versions.mix(GTF2TXDB.out.versions)

    emit:
        ch_bsgenome_file
        ch_bsgenome_name
        ch_txdb_file
        ch_versions
}
