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
        ch_sample_file
        ch_txdb
        ch_versions
    
    main:

        // CAGEr analysis steps
        if (params.bowtie2) {
            ch_data_type = Channel.of("bam")
        } else {
            ch_data_type = Channel.of("bigwig")
        }

        CAGER_READIN (
            ch_bsgenome_file,
            ch_bsgenome_name,
            ch_sample_file,
            ch_data_type
        )

        cager_rds = CAGER_READIN.out.rds
        ch_versions = ch_versions.mix(CAGER_READIN.out.versions)

        CAGER_TAG_QC(cager_rds, ch_txdb, ch_bsgenome_file, ch_bsgenome_name)
        ch_versions = ch_versions.mix(CAGER_TAG_QC.out.versions)

        CAGER_PREPROCESSING(cager_rds)
        clustered_cager_rds = CAGER_PREPROCESSING.out.rds
        ch_versions = ch_versions.mix(CAGER_PREPROCESSING.out.versions)

        CAGER_TAGCLUSTER_QC(clustered_cager_rds, ch_txdb, ch_bsgenome_file, ch_bsgenome_name)
        ch_versions = ch_versions.mix(CAGER_TAGCLUSTER_QC.out.versions)

        // TODO:
        // 1. consensus clusters
        // 2. track exports (what kinds?)
        // 3. expression profiling
        // 4. differential expression analysis
        // 5. shifting promoters
        // 6. enhancer calling

    emit:
        ch_versions

}
