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

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { CAGER } from '../modules/local/cager.nf'
include { DOWNLOAD_FASTA } from '../modules/local/downloadfasta.nf'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main.nf' 
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main.nf'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'

def multiqc_report = []

workflow CUSTOMCAGE {

    ch_versions = Channel.empty()

    if (!params.bsgenome) {
        exit 1, '--bsgenome (either a genome name from UCSC or a file path to a tar.gz archive) is not specified.'
    }

    if (params.input) {
        input_handler = file(params.input, checkIfExists: true)
    } else {
        exit 1, '--input (input samplesheet) is not specified.'
    }

    if (!params.fasta && !params.index) {
        bsgenome_name = file(params.bsgenome).name
        values = bsgenome_name.split('\\.')
        if (values[2] == "UCSC") {
            ucscid = values[3].split('_')[0]
            DOWNLOAD_FASTA( ucscid )
            ch_fasta = DOWNLOAD_FASTA.out.fasta
                .map{ it -> [[id: "FASTA"], it] }
        } else {
            exit 1, 'Reference fasta is not specified for a custom BSgenome.'
        }
    } else if (params.fasta && params.index) {
        exit 1, 'either --fasta or --index should be specified.'
    } else if (params.fasta) {
        fasta = [[[id: "FASTA"], params.fasta]]
        Channel
            .from( fasta )
            .map{ row -> [ row[0], file(row[1], checkIfExists: true) ] }
            .set{ ch_fasta }
    } else {
        index = [[[id: "INDEX"], params.index]]
        Channel
            .from( index )
            .map{ row -> [ row[0], file(row[1], checkIfExists: true) ] }
            .set{ ch_index }
    }

    INPUT_CHECK (
        input_handler
    )
    
    INPUT_CHECK.out.reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .map{ meta, fastq -> [ meta, fastq.flatten() ] }
        .set { ch_fastq }

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    CAT_FASTQ (
        ch_fastq
    ).reads.set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    TRIMGALORE (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    if (!params.index) {
        BOWTIE2_BUILD (
            ch_fasta
        )
        ch_index = BOWTIE2_BUILD.out.index
        //ch_index.view()
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())
    }

    ch_index1 = ch_index.map { it[1] }
    ch_index1.view()
    
    BOWTIE2_ALIGN (
        TRIMGALORE.out.reads,
        ch_index1, // ch_index
        false,
        false
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    SAMTOOLS_SORT (
        BOWTIE2_ALIGN.out.aligned
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // do not forget aligned reads must be processed with multiqc

    CAGER (
        params.bsgenome,
        SAMTOOLS_SORT.out.bam.collect()
    )
    ch_versions = ch_versions.mix(CAGER.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    workflow_summary    = WorkflowCustomcage.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCustomcage.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

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
