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

// Read deduplication parameters
params.dedup = false
params.dist = false

// HISAT2 parameters
params.hisat2 = false
params.gtf = "$projectDir/assets/NO_FILE_GTF"
params.splicesites = "$projectDir/assets/NO_FILE_SPLICESITES"
params.seq_center = false
params.save_unaligned = false
params.hisat2_build_memory = false

//TrimGalore! parameters
params.params_trimgalore = ''

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { CAGER } from '../modules/local/cager.nf'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main.nf' 
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main.nf'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as SORT_FOR_FIXMATE} from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP} from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_FIXMATE } from '../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_DEDUP } from '../modules/local/samtools_dedup.nf'
include { SAMTOOLS_STATS } from '../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_IDXSTATS } from '../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/samtools/flagstat/main.nf'
include { HISAT2_BUILD } from '../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN } from '../modules/nf-core/hisat2/align/main'

def multiqc_report = []

workflow CUSTOMCAGE {

    ch_versions = Channel.empty()
    ch_fasta = Channel.empty()

    if (!params.bsgenome) {
        exit 1, '--bsgenome (either a genome name from UCSC or a file path to a tar.gz archive) is not specified.'
    }

    if (params.input) {
        input_handler = file(params.input, checkIfExists: true)
    } else {
        exit 1, '--input (input samplesheet) is not specified.'
    }

    if (!params.fasta && !params.index) {
        exit 1, 'Reference FASTA file (--fasta) or genome index (--index) should be specified.'
    } else if (params.fasta && params.index) {
        exit 1, 'Only one of the two options, --fasta or --index, should be specified.'
    } else if (params.fasta) {
        Channel
            .fromPath(params.fasta)
            .set { ch_fasta }
    } else {
        Channel
            .fromPath(params.index)
            .set { ch_index }
    }

    if (params.dist) {
        if (!params.dedup) {
            exit 1, 'The --dist option can only be used with the --dedup option.'
        }
    }

    if (params.gtf != "$projectDir/assets/NO_FILE_GTF" & (!params.hisat2 || !params.fasta)) {
        exit 1, 'The --gtf option can only be used with the combination of the --hisat2 and --fasta options.'
    }

    if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" && !params.hisat2) {
        exit 1, 'The --splicesites option can only be used with the --hisat2 option.'
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

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    TRIMGALORE (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    if (params.hisat2) {            
        splice_sites_file = file(params.splicesites, checkIfExists: true)
        if (!params.index) {
            gtf_file = file(params.gtf, checkIfExists: true)
            HISAT2_BUILD (
                ch_fasta,
                [[id: "GTF"], gtf_file],
                [[id: "splice_sites"], splice_sites_file]
            )
            ch_index = HISAT2_BUILD.out.index
            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
        HISAT2_ALIGN (
            TRIMGALORE.out.reads,
            ch_index,
            splice_sites_file
        )
        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
        ch_aligned = HISAT2_ALIGN.out.bam
    } else {
        if (!params.index) {
            BOWTIE2_BUILD (
                ch_fasta
            )
            ch_index = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
        BOWTIE2_ALIGN (
            TRIMGALORE.out.reads,
            ch_index,
            false,
            false
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
        ch_aligned = BOWTIE2_ALIGN.out.aligned
    }

    if (params.dedup) {
        SORT_FOR_FIXMATE (
            ch_aligned
        )
        ch_versions = ch_versions.mix(SORT_FOR_FIXMATE.out.versions)

        SAMTOOLS_FIXMATE (
            SORT_FOR_FIXMATE.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions)
    }

    if (params.dedup) {
        ch_bam_to_sort = SAMTOOLS_FIXMATE.out.bam
    } else {
        ch_bam_to_sort = ch_aligned
    }

    SAMTOOLS_SORT (
        ch_bam_to_sort
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    if (params.dedup) {
        SAMTOOLS_DEDUP (
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_DEDUP.out.versions)

        SAMTOOLS_INDEX_DEDUP (
             SAMTOOLS_DEDUP.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions)
    }

    if (params.dedup) {
        ch_bam_bai = SAMTOOLS_DEDUP.out.bam.join(SAMTOOLS_INDEX_DEDUP.out.bai)
    } else {
        ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    }

    SAMTOOLS_STATS ( 
        ch_bam_bai, 
        ch_fasta.ifEmpty(file("$projectDir/assets/NO_FILE_FASTA", checkIfExists: true))
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    if (params.dedup) {
        ch_for_cager = SAMTOOLS_DEDUP.out.bam.collect()
    } else {
        ch_for_cager = SAMTOOLS_SORT.out.bam.collect()
    }

    CAGER (
        params.bsgenome,
        ch_for_cager
    )
    ch_versions = ch_versions.mix(CAGER.out.versions)

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
    if (params.hisat2) {
        ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})
    } else {
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})
    }
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})

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
