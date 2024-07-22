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

// input parameters
params.samplesheet = false

// TrimGalore! parameters
params.params_trimgalore = ''

// cutadapt parameters
params.nogtrim = false

// Read deduplication parameters
params.dedup = false
params.dist = false

// bowtie2 parameters
params.bowtie2 = false

// STAR parameters
params.gtf = "$projectDir/assets/NO_FILE_GTF"
params.chromsizes = "$projectDir/assets/NO_FILE_CHROMSIZES"
params.splicesites = "$projectDir/assets/NO_FILE_SPLICESITES"

// BSgenome parameters
params.bsgenome = false
params.forgeseed = false
params.sourcedir = false

include { INPUT_FROM_FOLDER } from '../subworkflows/local/input_from_folder.nf'
include { INPUT_FROM_SAMPLESHEET } from '../subworkflows/local/input_from_samplesheet.nf'
include { CAGER_BAM } from '../modules/local/cager_bam.nf'
include { CAGER_BIGWIG } from '../modules/local/cager_bigwig.nf'
include { FORGE_BSGENOME } from '../modules/local/forge_bsgenome.nf'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main.nf'
include { FASTQC } from '../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../modules/nf-core/multiqc/main.nf'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main.nf'
include { CUTADAPT } from '../modules/nf-core/cutadapt/main.nf'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main.nf' 
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main.nf'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main.nf' 
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main.nf'
include { SAMTOOLS_VIEW_MAPQ } from '../modules/nf-core/samtools/view_mapq/main.nf'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as SORT_FOR_FIXMATE} from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP} from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_FIXMATE } from '../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_DEDUP } from '../modules/local/samtools_dedup.nf'
include { SAMTOOLS_STATS } from '../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_IDXSTATS } from '../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/samtools/flagstat/main.nf'
include { UCSC_WIGTOBIGWIG } from '../modules/nf-core/ucsc/wigtobigwig/main.nf' 

def multiqc_report = []

workflow CUSTOMCAGE {

    ch_versions = Channel.empty()
    ch_fasta = Channel.empty()

    if (!params.bsgenome && (!params.forgeseed || !params.sourcedir)) {
        exit 1, 'Either the --bsgenome option or the following two options must be specified: --forgeseed, --sourcerdir.'
    } else if (params.bsgenome && (params.forgeseed || params.sourcedir)) {
        exit 1, 'The --bsgenome option and the following two options are mutually exclusive: --forgeseed, --sourcerdir.'
    }

    if (params.samplesheet) {
        input_handler = file(params.samplesheet, checkIfExists: true)
        INPUT_FROM_SAMPLESHEET (
            input_handler
        )

        INPUT_FROM_SAMPLESHEET.out.reads
            .map {
                meta, fastq ->
                    meta.id = meta.id.split('_')[0..-2].join('_')
                    [ meta, fastq ] }
            .groupTuple(by: [0])
            .map{ meta, fastq -> [ meta, fastq.flatten() ] }
            .set { ch_fastq }

        ch_versions = ch_versions.mix(INPUT_FROM_SAMPLESHEET.out.versions)
    } else if (params.infolder) {
        ch_fastq = INPUT_FROM_FOLDER(
            params.infolder
        )
    } else {
        exit 1, 'Provide input with the --samplesheet or the --infolder option.'
    }

    if (!params.fasta && !params.index) {
        exit 1, 'Reference FASTA file (--fasta) or genome index (--index) should be specified.'
    } else if (params.fasta && params.index) {
        exit 1, 'The --fasta and --index options are mutually exclusive.'
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
            exit 1, 'The --dist option requires the --dedup option.'
        }
    }

    if (params.gtf != "$projectDir/assets/NO_FILE_GTF" & !params.fasta) {
        exit 1, 'The --gtf option can only be used with the --fasta option.'
    }

    if (params.gtf != "$projectDir/assets/NO_FILE_GTF" & params.bowtie2) {
        exit 1, 'The --gtf option is mutually exclusive with the --bowtie2 option.'
    }

    if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" & !params.fasta) {
        exit 1, 'The --splicesites option can only be used with the --fasta option.'
    }

    if (params.splicesites != "$projectDir/assets/NO_FILE_SPLICESITES" && params.bowtie2) {
        exit 1, 'The --splicesites option is mutually exclusive with the --bowtie2 option.'
    }

    if (params.chromsizes != "$projectDir/assets/NO_FILE_CHROMSIZES" & params.bowtie2) {
        exit 1, 'The --chromsizes option is mutually exclusive with the --bowtie2 option.'
    }

    if (params.chromsizes == "$projectDir/assets/NO_FILE_CHROMSIZES" & !params.bowtie2) {
        exit 1, 'The use of the default mapper STAR requires the --chromsizes option.'
    }

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

    if (!params.nogtrim) {
        CUTADAPT (
            TRIMGALORE.out.reads
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
    }

    ch_reads_to_align = !params.nogtrim ? CUTADAPT.out.reads : TRIMGALORE.out.reads

    if (!params.bowtie2) {            
        if (!params.index) {
            gtf_file = file(params.gtf, checkIfExists: true)
            splice_sites_file = file(params.splicesites, checkIfExists: true)
            STAR_GENOMEGENERATE (
                ch_fasta,
                gtf_file,
                splice_sites_file
            )
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            
            ch_index = STAR_GENOMEGENERATE.out.index
        }

        STAR_ALIGN (
            ch_reads_to_align,
            ch_index
        )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        ch_aligned = STAR_ALIGN.out.bam_sorted

        ch_chrom_sizes = Channel.fromPath(params.chromsizes)

        UCSC_WIGTOBIGWIG (
            STAR_ALIGN.out.wigtobigwig,
            ch_chrom_sizes
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions)
    } else {
        if (!params.index) {
            BOWTIE2_BUILD (
                ch_fasta
            )
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
            
            ch_index = BOWTIE2_BUILD.out.index
        }

        BOWTIE2_ALIGN (
            ch_reads_to_align,
            ch_index,
            false,
            false
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

        SAMTOOLS_VIEW_MAPQ (
            BOWTIE2_ALIGN.out.aligned
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MAPQ.out.versions)

        ch_aligned = SAMTOOLS_VIEW_MAPQ.out.bam
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
            ch_bsgenome_file = file(params.bsgenome, checkIfExists: true)
            ch_bsgenome_name = ''
        } else {
            ch_bsgenome_file = file("$projectDir/assets/NO_FILE_BSGENOME")
            ch_bsgenome_name = params.bsgenome
        }
    } else {
        ch_bsgenome_file = FORGE_BSGENOME.out.bsgenome
        ch_bsgenome_name = ''
    }

    if (params.bowtie2) {
        if (params.dedup) {
            ch_for_cager = SAMTOOLS_DEDUP.out.bam.collect()
        } else {
            ch_for_cager = SAMTOOLS_SORT.out.bam.collect()
        }
        CAGER_BAM (
            ch_bsgenome_file,
            ch_bsgenome_name,
            ch_for_cager
        )
        ch_versions = ch_versions.mix(CAGER_BAM.out.versions)
    } else {
        ch_for_cager = UCSC_WIGTOBIGWIG.out.bw
            .map { it[1] }
            .collect()

        CAGER_BIGWIG (
            ch_bsgenome_file,
            ch_bsgenome_name,
            ch_for_cager
        )
        ch_versions = ch_versions.mix(CAGER_BIGWIG.out.versions)
    }

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
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
    if (!params.bowtie2) {
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})
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
