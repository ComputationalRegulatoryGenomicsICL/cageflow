/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog;
        paramsSummaryMap;
        validateParameters;
        paramsHelp;
        fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
// log.info logo + paramsSummaryLog(workflow) + citation

// WorkflowCustomcage.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
// include {  } from '../modules/nf-core/modules/ /main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CUSTOMCAGE {

    ch_versions = Channel.empty()
    // Channel.fromFilePairs("/Users/dbaranasic/data/playground/mock_fq/*_L00{1,2}_*.fastq.gz")
    // Channel
    //     // .fromFilePairs("/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/*_L00{1,2}_*.fastq.gz")
    //     // .fromPath("/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/*_L00{1,2}_*.fastq.gz")
    //     .fromFilePairs("/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/*_R{1,2}_*.fastq.gz")
    //     .set{ ch_reads_pe }
    // ch_reads_pe.view()

    // Channel
    //     .fromPath("/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/*_L00{1,2}*R1*.fastq.gz")
    //     .set{ ch_reads_se }
    // ch_reads_se.view()
        
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

    INPUT_CHECK (
        // ch_reads_pe
        // ch_reads
        ch_input
        // file(params.input)
    )

    // INPUT_CHECK.out.reads.view()
    
    INPUT_CHECK.out.reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    // ch_fastq.single.view()
    // ch_fastq.multiple.view()

    // ch_fastq.view()

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    ch_cat_fastq.view()
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    // if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

    // INPUT_CHECK (
    //     ch_input
    // )
    // .map {
    //     meta, fastq ->
    //         meta.id = meta.id.split('_')[0..-2].join('_')
    //         [ meta, fastq ] }
    // .groupTuple(by: [0])
    // .branch {
    //     meta, fastq ->
    //         single  : fastq.size() == 1
    //             return [ meta, fastq.flatten() ]
    //         multiple: fastq.size() > 1
    //             return [ meta, fastq.flatten() ]
    // }
    // .set { ch_fastq }


    // Create a new channel of metadata from a sample sheet
    // NB: `input` corresponds to `params.input` and associated sample sheet schema
    // ch_reads = Channel.fromSamplesheet(
    //     "input",
    //     parameters_schema: '/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/schema_input.json',
    //     skip_duplicate_check: false)
    // ch_reads = Channel.fromSamplesheet(
    //     "input")
    // ch_reads.view()

    //
    // MODULE: Run FastQC
    //

    // srr_paired = [
    // [
    // [
    //     id: "S10"
    // ],
    // "/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/S10_S6_L001_R1_001.fastq.gz",
    // "/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq/S10_S6_L001_R2_001.fastq.gz"
    // ],
    // ]

    // Channel
    //     .from( srr_paired )
    //     .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
    //     .set{ ch_srr_paired }

    // ch_srr_paired.view()

    // FASTQC (
    //     ch_reads_se
    //     // INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

/*

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCustomcage.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCustomcage.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
*/
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
