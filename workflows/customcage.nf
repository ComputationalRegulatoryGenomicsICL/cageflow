include { paramsSummaryLog;
          paramsSummaryMap;
          validateParameters;
          paramsHelp;
          fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// pipeline run settings
params.fullpipeline = true
params.maponly = false
params.cageronly = false

// input parameters
params.samplesheet = false

// TrimGalore! parameters
params.params_trimgalore = ''

// cutadapt parameters
params.nogtrim = false

// read deduplication parameters
params.dedup = false
params.dist = false

// genome annotation in GTF
params.gtf = "$projectDir/assets/NO_FILE_GTF"

// bowtie2 parameters
params.bowtie2 = false

// CAGEr markdown template location
params.markdown_path = "$projectDir/assets/cager_report.Rmd"

// BSgenome parameters
params.bsgenome = false
params.forgeseed = false
params.sourcedir = false

include { BIGWIG_INPUTS } from "../subworkflows/local/read_in_bigwigs.nf"
include { RELATIVISATION } from '../modules/local/make_paths_relative.nf'

include { PARAMETER_CHECKS } from '../subworkflows/local/parameter_checks.nf'
include { PREPROCESSING } from '../subworkflows/local/preprocessing.nf'
include { PREPARE_MAPPING_METADATA } from '../subworkflows/local/prepare_mapping_metadata.nf'
include { PREPARE_CAGER_METADATA } from '../subworkflows/local/prepare_cager_metadata.nf'
include { STAR_PROCESSING } from '../subworkflows/local/star_processing.nf'
include { BOWTIE2_PROCESSING } from '../subworkflows/local/bowtie2_processing.nf'
include { DEDUP } from '../subworkflows/local/deduplication.nf'
include { SAMTOOLS_PROCESSING } from '../subworkflows/local/samtools_processing.nf'
include { SUMMARY_STAT } from '../subworkflows/local/summary_statistics.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'
include { WRITE_SAMPLE_LIST } from '../modules/local/write_sample_list.nf'
include { CAGER } from '../subworkflows/local/cager_analysis.nf'

def multiqc_report = []

workflow CUSTOMCAGE {

    if (params.gtf) {
            ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
            // ch_pre_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
            // ch_gtf = sample_meta.combine(ch_pre_gtf)
            // ch_gtf = ch_genome_name.combine(ch_pre_gtf)
        } else {
            exit 1, "The --gtf argument is mandatory."
    }

    if (!params.maponly && !params.fullpipeline){
        if (!params.cager_sample_file ) {
            exit 1, 'Sample list file is mandatory if mapping is not done within the pipeline.'
        }

        ch_cager_sample_file = Channel.fromPath(params.cager_sample_file)
        bigwig_files_ch = BIGWIG_INPUTS(ch_cager_sample_file).collect()
        merged_sample_file = RELATIVISATION(ch_cager_sample_file)

    }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (params.maponly || params.fullpipeline) {

        ch_fasta = Channel.empty()
        ch_index = Channel.empty()

        PARAMETER_CHECKS(ch_fasta, ch_index, ch_versions)

        ch_fasta = PARAMETER_CHECKS.out.ch_fasta
        ch_index = PARAMETER_CHECKS.out.ch_index
        ch_fastq = PARAMETER_CHECKS.out.ch_fastq
        ch_versions = PARAMETER_CHECKS.out.ch_versions

        PREPROCESSING(ch_fastq, ch_versions, ch_multiqc_files)

        ch_reads_to_align = PREPROCESSING.out.ch_reads_to_align
        ch_multiqc_files = PREPROCESSING.out.ch_multiqc_files
        ch_versions = PREPROCESSING.out.ch_versions
        
        PREPARE_MAPPING_METADATA( ch_fasta, ch_versions )
        ch_chrom_sizes = PREPARE_MAPPING_METADATA.out.ch_chrom_sizes
        ch_fasta = PREPARE_MAPPING_METADATA.out.ch_fasta
        ch_versions = PREPARE_MAPPING_METADATA.out.ch_versions

        if (params.bowtie2) {            
            BOWTIE2_PROCESSING(ch_reads_to_align, ch_fasta, ch_index, ch_multiqc_files, ch_versions)
            
            ch_aligned = BOWTIE2_PROCESSING.out.ch_aligned
            ch_multiqc_files = BOWTIE2_PROCESSING.out.ch_multiqc_files
            ch_versions = BOWTIE2_PROCESSING.out.ch_versions
            // NOTE: placeholder so that the channel is not empty
            // it will be replaced in SAMTOOLS_PROCESSING
            ch_for_cager = ch_aligned

        } else {
            STAR_PROCESSING(ch_reads_to_align, ch_fasta, ch_index, ch_gtf, ch_chrom_sizes, ch_multiqc_files, ch_versions)

            ch_for_cager = STAR_PROCESSING.out.bigwig_ch_for_cager
            ch_aligned = STAR_PROCESSING.out.ch_aligned
            ch_multiqc_files = STAR_PROCESSING.out.ch_multiqc_files
            ch_versions = STAR_PROCESSING.out.ch_versions
        }

        if (params.dedup) {
            DEDUP(ch_aligned, ch_versions, ch_for_cager)

            ch_for_cager = DEDUP.out.ch_for_cager
            ch_bam_bai = DEDUP.out.ch_bam_bai
            ch_versions = DEDUP.out.ch_versions
        } else {
            SAMTOOLS_PROCESSING(ch_aligned, ch_versions, ch_for_cager)

            ch_for_cager = SAMTOOLS_PROCESSING.out.ch_for_cager
            ch_bam_bai = SAMTOOLS_PROCESSING.out.ch_bam_bai
            ch_versions = SAMTOOLS_PROCESSING.out.ch_versions
        }

        SUMMARY_STAT(ch_bam_bai, ch_fasta, ch_multiqc_files, ch_versions)

        ch_multiqc_files = SUMMARY_STAT.out.ch_multiqc_files
        ch_versions = SUMMARY_STAT.out.ch_versions

        bigwig_files_ch = ch_for_cager.map{ meta, paths ->
            file1 =  paths[0]
            file2 = paths[1]
            [file1, file2]}
            .collect()

        // add new_name equal sample ID
        ch_sample_files = WRITE_SAMPLE_LIST(ch_for_cager)
        // add new_name
        def header = "id,single_end,path"

        ch_collected = ch_sample_files
        .reduce( header ) { acc, table_line ->
            acc + '\n' + table_line.readLines()[0]}

        // sorting samples alphabetically
        merged_sample_file = ch_collected.collectFile(
            name: "sample_list.csv",
            newLine: true,
            sort: { file -> file.text })

        CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        workflow_summary    = WorkflowCustomcage.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)
        methods_description    = WorkflowCustomcage.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
    
    if (params.cageronly || params.fullpipeline) {

        PREPARE_CAGER_METADATA( ch_gtf, ch_versions )
        ch_bsgenome_file = PREPARE_CAGER_METADATA.out.ch_bsgenome_file
        ch_bsgenome_name = PREPARE_CAGER_METADATA.out.ch_bsgenome_name
        ch_txdb_file = PREPARE_CAGER_METADATA.out.ch_txdb_file
        ch_versions = PREPARE_CAGER_METADATA.out.ch_versions

        CAGER(
            ch_bsgenome_file,
            ch_bsgenome_name,
            merged_sample_file,
            bigwig_files_ch,
            ch_txdb_file,
            ch_versions
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.dump_parameters(workflow, params)
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
