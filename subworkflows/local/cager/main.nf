//
// CAGEr analysis steps
//
include { CAGER_READIN } from '../../../modules/local/cager/readin/main.nf'
include { CAGER_TAG_QC } from '../../../modules/local/cager/tag_qc/main.nf'
include { CAGER_PROCESSING } from '../../../modules/local/cager/processing/main.nf'
include { CAGER_TAGCLUSTER_QC } from '../../../modules/local/cager/tagcluster_qc/main.nf'
include { CAGER_REPORT } from "../../../modules/local/cager/report/main.nf"
include { CAGEFIGHTR_ENHANCER_CALLING } from '../../../modules/local/cagefightr/enhancer_calling/main.nf'


workflow CAGER {

    take:
        ch_bsgenome_file
        ch_bsgenome_name
        ch_sample_file
        ch_collected
        ch_txdb
        ch_versions

    main:

        // CAGEr analysis steps
        if (params.bowtie2) {
            ch_data_type = Channel.of("bam")
        } else {
            ch_data_type = Channel.of("bigwig")
        }

        sample_table = ch_sample_file
            .splitCsv( header:true , sep:',')
            .map { create_mapping_channel(it) }
            .collect()

        CAGER_READIN (
            ch_bsgenome_file,
            ch_bsgenome_name,
            sample_table,
            ch_data_type,
            ch_collected
        )

        cager_rds = CAGER_READIN.out.rds
        ch_versions = ch_versions.mix(CAGER_READIN.out.versions)

        CAGER_TAG_QC(cager_rds, ch_txdb, ch_bsgenome_file, ch_bsgenome_name)
        annotated_cager_rds = CAGER_TAG_QC.out.cager_rds
        ch_versions = ch_versions.mix(CAGER_TAG_QC.out.versions)
        // tag region annotation, correlations heatmap
        tra_ch_tss = CAGER_TAG_QC.out.plots
        tag_corr_data = CAGER_TAG_QC.out.correlation_rds

        CAGER_PROCESSING(annotated_cager_rds, ch_bsgenome_file, ch_bsgenome_name, ch_txdb)
        clustered_cager_rds = CAGER_PROCESSING.out.rds
        ch_versions = ch_versions.mix(CAGER_PROCESSING.out.versions)
        // reverse cumulative, iterquartile width, ctss counts
        ch_preproc_res = CAGER_PROCESSING.out.results

        CAGER_TAGCLUSTER_QC(clustered_cager_rds, ch_txdb, ch_bsgenome_file, ch_bsgenome_name)
        ch_versions = ch_versions.mix(CAGER_TAGCLUSTER_QC.out.versions)
        // tagcluster annotations, nucleotide frequencies, dinucleotide frequencies, TSSlogo
        ch_tagc_plots = CAGER_TAGCLUSTER_QC.out.plots
        ch_plots = ch_tagc_plots
        tc_corr_data = CAGER_TAGCLUSTER_QC.out.correlation_rds

        // enhancer calling
        CAGEFIGHTR_ENHANCER_CALLING(
            clustered_cager_rds,
            ch_txdb)
        ch_versions = ch_versions.mix(CAGEFIGHTR_ENHANCER_CALLING.out.versions)
        // enhancer calling plots
        enhancer_plots = CAGEFIGHTR_ENHANCER_CALLING.out.plots

        ch_template = Channel.fromPath(params.markdown_path)

        ch_html = CAGER_REPORT(
            ch_template,
            tra_ch_tss,
            tag_corr_data,
            ch_preproc_res,
            ch_tagc_plots,
            tc_corr_data,
            enhancer_plots,
            cager_rds)

        ch_plots = ch_plots.mix(enhancer_plots)

    emit:
        ch_plots
        ch_html
        ch_versions

}

def create_mapping_channel(LinkedHashMap row) {
    id = row.id
    single_end = row.single_end
    str1_bw = row.path.split(" ")[0].minus('[')
    str2_bw = row.path.split(" ")[1].minus(']')
    new_name = row.new_name

    return [id, single_end, str1_bw, str2_bw, new_name]
}
