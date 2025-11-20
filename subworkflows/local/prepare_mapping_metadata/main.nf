//
// Subworkflow to get the chromsizes
//

include { SAMTOOLS_FAIXD } from '../../../modules/nf-core/samtools/faidx/main.nf'

workflow PREPARE_MAPPING_METADATA {

    take:
        ch_fasta
        ch_versions

    main:

        // prepare chromosome sizes
        if (params.fasta) {

            chrom_size_fa = ch_fasta.map{ meta, fasta ->
                def new_meta = [:]
                new_meta.id = "sizes"
                fasta = fasta
                [new_meta, fasta, true]
            }.unique()
            SAMTOOLS_FAIXD( chrom_size_fa )
            ch_chrom_sizes = SAMTOOLS_FAIXD.out.sizes

            ch_versions = ch_versions.mix(SAMTOOLS_FAIXD.out.versions)
        } else { // a genome index was provided instead
            ch_chrom_sizes = Channel.of([
                [id:"sizes"],
                [file( params.index + '/chrNameLength.txt' )]
            ])
        }

    emit:
        ch_chrom_sizes
        ch_fasta
        ch_versions
}
