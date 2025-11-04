# ComputationalRegulatoryGenomicsICL/customcage: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

**ComputationalRegulatoryGenomicsICL/customcageq** is a Nextflow pipeline to process CAGE sequencing data from raw reads to the identification of consensus clusters and enhancers.
The pipeline is specifically designed to be used for assessing the quality of the CAGE data, and calling promoters and enhancers with good default values that can be later fine tuned as needed.
The pipeline also supports the re-analysis of the reads with updated parameters.


### Input

#### Parameter file

Customizing the analysis the pipeline requires about 40 different parameters. These can be defined in a single `params.yaml` file.
An example is shown below.

```
# pipeline parameters
fullpipeline: true
maponly: false
cageronly: false
gtf: "testdata/sacCer3_genome/sacCer3.ensGene.gtf"

# preprocessing parameters
samplesheet: "docs/examples/samplesheet_sacer_pe.csv"
infolder:
sample_name_fields:
# mapping parameters
genome_name: "sacCer3"
fasta: "testdata/sacCer3_genome/sacCer3.fa"
index: "testdata/sacCer3_genome/sacCer3_star_index/"
seq_platform: "illumina"
seq_center: false
# whether to take the uniquely mapped only
unique_only: true
# whether to remove reads that do not start with G
remove_non_g: false 

# CAGEr parameters
cager_sample_file: "docs/examples/sample_list.csv" # with sorted list of bigwigs, if mapping is run elsewhere
# BSgenome
forgeseed:
sourcedir:
bsgenome: "BSgenome.Scerevisiae.UCSC.sacCer3"
# parameter for correlation calculation
corrplot_tagCountThreshold: 1
# parameters for normalization
norm_range_min: 5
norm_range_max: 10000
norm_method: "powerLaw"
alpha:
T_norm: 1000000
# parameters for tag clustering
sample_num_thr: 1
ctss_thr: 1
distclu_maxDist: 20
keepSingletonsAbove: 5
iq_low: 0.1
iq_high: 0.9
# plotting for tagclusters QC
iqw_tpm_threshold: 3
tssregion_up: -3000
tssregion_down: 3000
tsslogo_upstream: 35
# parameters for consensus clusters
consensus_thr: 2
consensus_dist: 100
# parameters for enhancer calling
cfBalanceThreshold: 0.95
unexpressed: 0
minSamples: 0
```

where the pipeline parameters that should be provided in all runs are
* `fullpipeline`, `maponly`, and `cageronly` define the mode of the pipeline run. With the `fullpipeline` option, every step will run. `maponly` stops at creating `bigwig` and/or `bam` files, as well as a MultiQC report. `cageronly` starts from `bigwig` or `bam` files and finishes with a CAGEr report.
* `gtf` specifies a GTF file with a whole-genome annotation to create a STAR index (if a FASTA file is provided) and/or to automatically create a TxDb annotation object to annotate tag clusters.

The parameters specific to mapping, can be left empty when running in `cageronly` mode:
* `samplesheet` specifies the input CSV samplesheet. This option is mutually exclusive with `infolder`.
* `infolder` specifies the input directory with FASTQ files (stored together for all samples or located in per-sample subdirectories). This option is mutually exclusive with `samplesheet`, and may be used together with `sample_name_fields`.
* `sample_name_fields` is a supporting parameter for `infolder` in case your sample name has underscore(s) in it. By default, only the first part of the string before the first underscore is taken for samplename. If you have more, like `my_sample_name_S1_L001_R1_001.fastq.gz`, with this parameter you may specify *how many underscore separated fields* the sample name has in the filename. In the `my_sample_name_S1_L001_R1_001.fastq.gz` example, this parameter should be = 3.
* `genome_name` specifies the name of the reference genome. It is used as meta information
* `fasta` specifies a FASTA file containing a reference genome. This option is mandatory, unless `index` is set.
* `index` specifies a directory with a genome index (`bowtie2` or `STAR`). This is a mandatory option, unless `fasta` is set.
* `seq_platform` specifies the sequencing platform used. Required for mapping with `STAR`.
* `seq_center` specifies the name of the sequencing center. Required for mapping with `STAR`.
* `unique_only` specifies if only uniquely mapped reads are considered for downstream analysis. Required for mapping with `STAR`. Not considered when using `bowtie2`.
* `remove_non_g` specifies whether to keep only those reads that start with `G` base, as expected after the CAGE protocol. This step is expected to remove about 15-20% of the reads that would likely be non canonical initiators, but they might take part in ohter biological processes.


The parameters specific to CAGEr and CAGEfightR analysis, can be left empty when running in `maponly` mode:
* `cager_sample_file` specifies the input CSV samplesheet including the name of the samples, their pairedness status, and the location of bigwigs. Optionally, a fourth column is used that sepcifies which samples should be removed, merged, or kept as is. This is achieved by checking if the row is empty (remove), its content is unique (keep as is), or shared with another sample (merge).
* `datatype` specifies the type of mapped input in the CSV samplesheet: `bam` or `bigwig` is accepted. Default is `bigwig`.
* `forgeseed` specifies a seed file for BSgenome forging (see the [Advanced BSgenomeForge usage vignette](https://bioconductor.org/packages/release/bioc/vignettes/BSgenomeForge/inst/doc/AdvancedBSgenomeForge.pdf) for details). The seed file should not contain the `seqs_srcdir` field (instead, the absolute or relative path to the source directory is set with the `sourcedir` option, see below). This option requires `sourcedir` and is mutually exclusive with `bsgenome`.
* `sourcedir` specifies a directory containing either a set of FASTA files, one per reference chromosome, or a 2bit file for the whole reference genome. See the [Advanced BSgenomeForge usage vignette](https://bioconductor.org/packages/release/bioc/vignettes/BSgenomeForge/inst/doc/AdvancedBSgenomeForge.pdf) for details. The seed file should be written according to the contents of this directory. This option requires `forgeseed` and is mutually exclusive with `bsgenome`.
* `bsgenome` specifies the BSgenome R package to use. If it is a file name (which should have a full path and the `.tar.gz` extension), then the package will be taken from the specified location; otherwise, the pipeline will try to install a BSgenome R package with the name `bsgenome.package` on the fly (see examples below). This option is mutually exclusive with `forgeseed` and `sourcedir`.
* `corrplot_tagCountThreshold` is a threshold above which (raw and normalized) CTSS are considered for the correlation plot.
* `norm_method` is the method used for normalizing the samples. Options are `simpleTpm` to covert tag counts to tags per million, `powerLaw` to normalize to a reference power-law distribution, or `none` to keep using the raw tag counts in downstream analyses. Case sensitive.
* `norm_range_min` and `norm_range_max` defines the lower and upper thresold for fitting the power-law distribution and calculate the slope for normalization. Only used when `norm_method` is `powerLaw`.
* `alpha` user specified alpha, the `-1 *` fitted slope in the log-log representation of the power-law distribution. If none, the average across samples is calculated and used. Considered for `powerLaw` normalization only. Note: this value cannot be less than 1.05, if the user provides a lower value the code automatically assumes 1.05.
* `T_norm` total number of CAGE tags in the reference power-law distribution. Setting `T = 10^6` results in normalized values that correspond to tags per million. Considered for `powerLaw` normalization only.
* `sample_num_thr` and `ctss_thr` are parameters for filtering low expressed CTSS before clustering. `ctss_thr` specifies the lower threshold above which CTSS are considered, and `sample_num_thr` specifies the number of samples where this threshold should be passed.  
* `distclu_maxDist` specifies the maximum distance for distance-based clustering (distclu).
* `keepSingletonsAbove` defines the tpm threshold above which even a single CTSS is kept during clustering.
* The `iq_low` and `iq_high` parameters are used to define the lower and upper quantile boundaries of the interquartile range within which the majority of the signal lies.
* `iqw_tpm_threshold` is a threshold above which tag clusters are considered for the interquartile width distribution plot.
* `tssregion_up` and `tssregion_down` are used for annotation with ChIPseeker. These correspond to the upstream and downstream distance to consider into TSS region for ChIPseeker annotation.
tssregion_up should be negative, tssregion_down should be positive.
* `tsslogo_upstream` is used for plotting the TSS logos. This parameter specifies the number of bases to inlcude upstream of the TSS.
* `consensus_thr` and `consensus_dist` are used for defining the consensus clusters. `consensus_thr` specifies the TPM threshold above which tag clusters are considered for consensus clusters, and  `consensus_dist` define the maximum distance between the interquartile ranges of tag clusters to be joined together into consensus clusters.
* `cfBalanceThreshold` is used for enhancer calling. It defines the balance threshold above which bidirectionality is considered balanced.
* `unexpressed` and `minSamples` are used for selecting only supported enhancers. `unexpressed` is a non inclusive lower TPM boundary for expression when calculating support of enhancers. `minSamples` is a non-inclusive lower boundary for the number of samples where the clusters should show bidirectionality.


* Additional arguments can be:
    * `params-trimgalore 'params'` specifies any options that can be passed to `TrimGalore!`. This option is useful for any non-standard read processing (for example, for CAGEscan reads that require the removal of a fixed number of nucleotides from the 5'-ends of the forward and reverse reads ([Bertin et al., 2017](https://www.nature.com/articles/sdata2017147))). The string with the parameters for `TrimGalore!` must be surrounded by single quotes.
    * `nogtrim` makes the pipeline skip the G-trimming step. This option is useful for processing non-CAGE data (for example, CAGEscan reads which do not seem to require trimming of a 5'-`G` ([Bertin et al., 2017](https://www.nature.com/articles/sdata2017147))). This option can be used together with `params-trimgalore` (see an example below).
    * `bowtie2` switches the aligner from `STAR` to `bowtie2`. This option is compatible with either `index` or `fasta`.
    * `dedup` switches on PCR duplicate removal (not shown in the pipeline map above and is switched off by default).
    * `dist L` sets an optical duplicate distance `L` to remove optical duplicates, in addition to PCR duplicates (see [`samtools markdup`](https://www.htslib.org/doc/samtools-markdup.html), option `-d`). This option requires `dedup`.
    * `-w` is a Nextflow option that specifies a path to the Nextflow work directory.

*Note: All pipeline options may be provided to the nextflow command starting with a double dash (`--`). All Nextflow options start with a single dash (`-`).*

#### Complete pipeline

To run the complete pipeline, starting from raw reads, the input is either single-end (SE) or paired-end (PE) raw CAGE reads. **Only one type of reads (either SE or PE) can be used in one run of the pipeline.** The user can list read files in a samplesheet or provide a path to a directory containing the files (stored all together or in subdirectories within the provided folder).

##### Samplesheet
The samplesheet has to be a comma-separated file with 3 columns and a header row. A PE example is shown below (this file and a SE version can be found at `docs/examples/samplesheet_sacer_\[pe/se\].csv`). It is recommended to use absolute paths for the input files. 
If this option is used, the input is defined with the `--samplesheet` parameter.

```
sample,fastq_1,fastq_2,single_end
S10,testdata/sacCer_fastq/pe/S10_S6_L001_R1_001.fastq.gz,testdata/sacCer_fastq/pe/S10_S6_L001_R2_001.fastq.gz,False
S14,testdata/sacCer_fastq/pe/S14_S8_L001_R1_001.fastq.gz,testdata/sacCer_fastq/pe/S14_S8_L001_R2_001.fastq.gz,False
S14,testdata/sacCer_fastq/pe/S14_S8_L002_R1_001.fastq.gz,testdata/sacCer_fastq/pe/S14_S8_L002_R2_001.fastq.gz,False
```

where
* `sample` is a unique identifier of a sample;
* `fastq_1` (and `fastq_2` in the case of paired-end reads) is a full path to the read libraries. In case of paired-end reads, `fastq_1` contains the full path to forward reads, while `fastq_2` contains the full path to reverse reads. One sample can be represented by more than one library if lanes are stored separately;
* `single_end` should be set to `True` for single-end reads and to `False` for paired-end reads.

For paired-end reads, `fastq_2` should contain the full path to reverse reads, while `single_end` should be set to `False`.

- [Samplesheet](samplesheet.md)
  - Additional information on the samplesheet.


##### Starting from a folder path

If a foldername including fastq files is provided, the `--infolder` parameter should be selected.
In this case, the software assumes that the file name follows Illumina naming conventions and has the following structure: `<sample name>_<sample number>_L00<lane number>_<R1/R2/1/2>_001.fastq.gz`.
From the existence of `R2` (or `2`) values in the expected position, the software assigns the input to be paired end or single end. 
By default, it takes the first value before the underscore which should be the *sample name*.
However, it can be overwritten in case the *sample name* itself contains underscore values.
The additional parameter `sample_name_fields` should be set to how many underscore separated parts the *sample name* has.
*Note, that if the sample name had any dash (`-`) it is converted to underscore (`_`) to ensure compatibility with subsequent processes.*

#### CAGEr subpipeline

To run the analysis with CAGEr with already existing bigwig or bam files, the input is another sample sheet with 4 columns and a header row. 
the 4th column, `new_name` defines which samples should be merged (provide the same names for the samples), dropped (leave the field empty), or kept as is (match the with the id field).
An example is shown below (this file can be found at `docs/examples/samplesheet_cager.csv`). It is recommended to use absolute paths for the input files.

```
id,single_end,path,new_name
S10,false,[testdata/bigwig/S10.Signal.UniqueMultiple.str1.out.wig.bw testdata/bigwig/S10.Signal.UniqueMultiple.str2.out.wig.bw],S10
S14,false,[testdata/bigwig/S14.Signal.UniqueMultiple.str1.out.wig.bw testdata/bigwig/S14.Signal.UniqueMultiple.str2.out.wig.bw],S14
```

This file may be created in two steps as shown below.

```
ls /path/to/bigwig_folder/*bw > bigwigs.txt`.
./bin/sample_list_from_bigwig_list.py --filepath bigwigs.txt --singleend ["true" or "false"]
```

First, a simple txt file listing all bigwigs is created. Next, the custom script is run specifying the conversion which will create a file called `sample_list.csv`.
This script accepts two additional parameters: `--delimiter` and `--field` specifying if any part of the filename should be excluded. This might be necessary because a filename should not start with a number, as filename will be propagated as sample name and R is sensitive to names starting with a number.

## Running the pipeline

After you cloned this repository, the typical command for running the pipeline is as follows:

```bash
nextflow run customcageq/main.nf -params-file params.yaml -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir, default is results/)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Additional information

You can change the maximum number of instances of the same process that can run in parallel (by default, the maximum number of instances of the same process equals 2).
To do this, change the value of the `maxForks` parameter in `conf/base.config`.
<!-- TODO: check if this can be overwritten with the -c option in nextflow, that would be cleaner -->
Limiting this number makes sure that the pipeline does not try to obtain all available resources.

By default, up to 10 different mapping loci are allowed for each read but only uniquely mapping reads are selected to create bigWig files.
If you would like to work with multimappers (or change mapping parameters in any other way), please amend STAR mapping options in `conf/modules.config` and/or replace the default `true` value for the `unique_only` pipeline option with `false` (in which case bigWig files will be created using the full set of alignments, not only unique ones).
However, if multimapping reads are allowed, then bigWig files will contain non-integer counts for positions where multimappers align.
This is due to the fact that STAR splits the count of 1 between all alignments of the same read.

- [nextflow_usage](docs/nextflow_usage.md)
  - Nextflow specific information


## Examples

### Paired-end reads with locally stored STAR index and BSgenome

Call TSSs from the test yeast paired-end CAGE reads using the locally stored test STAR index and the `BSgenome.Scerevisiae.UCSC.sacCer3` R package. The package will be automatically downloaded and installed within the CAGEr container on the fly and will be used there with CAGEr. To run this example, the user needs to provide full paths to the test FASTQ files in `samplesheet_sacer_pe_template.csv` and the path to a "scratch" (or any other convenient) storage space for the Nextflow work directory:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer3 \
    --gtf customcageq/assets/sacCer3_genome/sacCer3.ensGene.gtf \
    --index customcageq/assets/sacCer3_genome/sacCer3_star_index/ \
    --samplesheet customcageq/assets/samplesheet_sacer_pe_template.csv \
    -profile singularity \
    -w /path/to/scratch/work
```

This example represents a typical use case for processing CAGE data from an organism with an available, locally stored, STAR genome index and a corresponding BSgenome package available in Bioconductor. For example, this is a use case for human CAGE data processing with the hg38 or T2T-CHM13 assembly. Instead of the `singularity` profile, one may use their institution's Nextflow profile (see [publicly available institutional Nextflow profiles](https://nf-co.re/configs)).

### Single-end reads with locally stored FASTA and GTF files

Call TSSs from the test yeast single-end CAGE reads (stored in per-sample subdirectories of an input directory) using locally stored FASTA and GTF files (and an optional file with splice junctions) for `STAR` index generation on the fly, as well as a locally stored seed file and a source directory to build a BSgenome R package. The package will be automatically installed within the CAGEr container from the autogenerated `.tar.gz` archive and used with CAGEr. To run this example, the user needs to provide full paths to the test FASTQ files in `samplesheet_sacer_se_template.csv`, as well as paths to a locally stored seed file and a source directory:

```bash
nextflow run customcageq/main.nf \
    --forgeseed /path/to/bsgenome_forging.seed \
    --sourcedir /path/to/seqs_srcdir \
    --gtf customcageq/assets/sacCer3_genome/sacCer3.ensGene.gtf \
    --fasta customcageq/assets/sacCer3_genome/sacCer3.fa \
    --gtf customcageq/assets/sacCer3_genome/sacCer3.ensGene.gtf \
    --infolder customcageq/assets/sacCer_fastq/pe_per_sample/ \
    -profile singularity
```

where the `pe_per_sample` input directory has the following structure:

```bash
customcageq/assets/sacCer_fastq/pe_per_sample/
├── S10
│   ├── S10_S6_L001_R1_001.fastq.gz
│   └── S10_S6_L001_R2_001.fastq.gz
└── S14
    ├── S14_S8_L001_R1_001.fastq.gz
    ├── S14_S8_L001_R2_001.fastq.gz
    ├── S14_S8_L002_R1_001.fastq.gz
    └── S14_S8_L002_R2_001.fastq.gz
```

This example may suit for the processing of CAGE data from a new species for which the user has to build ("forge") a BSgenome package by themselves. After forging the BSgenome once, the user can copy the resulting `.tar.gz` file from `results/bsgenome` and reuse it in subsequent runs of the pipeline by setting the `--bsgenome` option, for example: `--bsgenome /path/to/bsgenome/BSgenome.Scerevisiae.UCSC.sacCer3_1.4.0.tar.gz`.

### CAGEscan libraries

Call TSSs from FANTOM5 CAGEscan libraries (see, for example, [CAGEscan datasets from human primary cells by FANTOM5](https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.primary_cell.CAGEScan/)). These libraries require trimming of 9 nt from the 5'-ends of the forward reads and of 6 nt from the 5'-ends of the reverse reads and do not seem to require separate G-trimming ([Bertin et al., 2017](https://doi.org/10.1038/sdata.2017.147)). To run this example, the user needs to generate the `fantom5_cagescan_pe.csv` input table (see above) and provide a path to it and to the `STAR` index of the T2T-CHM13 v2.0 human genome assembly:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 \
    --gtf /path/to/chm13_t2t.gtf \
    --index /path/to/chm13_t2t_v2.0_star_index \
    --params-trimgalore '--clip_R1 9 --clip_R2 6' \
    --nogtrim \
    --samplesheet /path/to/fantom5_cagescan_pe.csv \
    -profile singularity
```

### **Not recommended.** Paired-end reads with bowtie2 and duplicate removal

Call TSSs from the test yeast paired-end CAGE reads using `bowtie2` for read mapping. Additionally, remove PCR duplicates and optical duplicates at a maximum distance 100 (see [`samtools markdup`](https://www.htslib.org/doc/samtools-markdup.html)) before doing alignment QC and TSS calling:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer3 \
    --gtf customcageq/assets/sacCer3_genome/sacCer3.ensGene.gtf \
    --bowtie2 \
    --index customcageq/assets/sacCer3_genome/sacCer3_bowtie2_index \
    --dedup \
    --dist 100 \
    --samplesheet customcageq/assets/samplesheet_sacer_pe_template.csv \
    -profile singularity
```

This example collects options that are **not recommended** but retained just in case. Using `bowtie2` does not allow accounting for splicing and makes CAGEr use BAM files, which slows the creation of the CAGEexp object considerably. Also, read deduplication is not recommended because CAGE reads, by design, come only from transcripts, with one read coming from the 5'-end of the transcript, which increases the probability of true duplicates, in comparison to whole-genome libraries, like ChIP-seq or ATAC-seq. The `--bowtie2` option can also be used with the `--fasta` option to build a `bowtie2` genome index on the fly, while the `--dedup` and `--dist` options can be used with the default `STAR` mapping.


## Credits

**ComputationalRegulatoryGenomicsICL/customcageq** has been developed by Sviatoslav Sidorov ([@sidorov-si](https://github.com/sidorov-si)), Katalin Ferenc ([@ferenckata](https://github.com/ferenckata)), Damir Baranasic ([@da-bar](https://github.com/da-bar)), Elena Gómez-Marín ([@ElenaGoMa](https://github.com/ElenaGoMa)), and Pavel Nikitin ([@nikitin-p](https://github.com/nikitin-p)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  ComputationalRegulatoryGenomicsICL/customcage for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

