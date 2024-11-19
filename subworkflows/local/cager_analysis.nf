// 
// CAGEr analysis steps
// 
include { CAGER_READIN } from '../../modules/local/cager_readin.nf'
include { CAGER_TAG_QC } from '../../modules/local/cager_tag_qc.nf'
include { CAGER_PREPROCESSING } from '../../modules/local/cager_preprocessing.nf'
include { CAGER_TAGCLUSTER_QC } from '../../modules/local/cager_tagcluster_qc.nf'


workflow CAGER {

    take:
        ch_bsgenome_file
        ch_bsgenome_name
        ch_for_cager
        bigwig_ch_for_cager
        ch_versions
    
    main:

        // CAGEr analysis steps
        if (params.bowtie2) {
            ch_data_type = Channel.of("bam")
            ch_data_in = ch_for_cager
        } else {
            ch_data_type = Channel.of("bigwig")
            ch_data_in = bigwig_ch_for_cager
        }

        ch_sample_file = WRITE_SAMPLE_LIST(ch_data_in)

        CAGER_READIN (
            ch_bsgenome_file,
            ch_bsgenome_name,
            ch_sample_file,
            ch_data_type
        )

        cager_rds = CAGER_READIN.out.rds
        ch_versions = ch_versions.mix(CAGER_READIN.out.versions)


        // CAGER_TAG_QC(cager_rds, ch_txdb)
        // ch_versions = ch_versions.mix(CAGER_TAG_QC.out.versions)

        // CAGER_PREPROCESSING(cager_rds)
        // clustered_cager_rds = CAGER_PREPROCESSING.out.rds
        // ch_versions = ch_versions.mix(CAGER_PREPROCESSING.out.versions)

        // CAGER_TAGCLUSTER_QC(clustered_cager_rds)
        // ch_versions = ch_versions.mix(CAGER_TAGCLUSTER_QC.out.versions)

    emit:
        ch_versions

}
