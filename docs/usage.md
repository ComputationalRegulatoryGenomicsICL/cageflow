# ComputationalRegulatoryGenomicsICL/customcage: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

**ComputationalRegulatoryGenomicsICL/customcageq** is a Nextflow pipeline to process CAGE sequencing data from raw reads to the identification of consensus clusters and enhancers. The pipeline is specifically designed to be used for assessing the quality of the CAGE data, and calling promoters and enhancers with good default values that can be later fine tuned as needed. The pipeline also supports the re-analysis of the reads with updated parameters.


### Input

#### Parameter file

To customize the analysis the pipeline takes about 40 different parameters, which can be defined in a single `params.yaml` file.
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
star_ignore_sjdbgtf: false
seq_platform: "illumina"
seq_center: false
# whether to take the uniquely mapped only
unique_only: false
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
```

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

- [Samplesheet](samplesheet.md)
  - Additional information on the samplesheet.


##### Starting from a folder path
If a foldername including fastq files is provided, the `--infolder` parameter should be selected.
In this case, the software assumes that the file name follows Illumina naming conventions and has the following structure: `<sample name>_<sample number>_L00<lane number>_<R1/R2/1/2>_001.fastq.gz`.
From the existence of `R2` (or `2`) values in the expected position, the software assigns the input to be paired end or single end. 
By default, it takes the first value before the underscore which should be the *sample name*.
However, it can be overwritten in case the *sample name* itself contains underscore values.
The additional parameter `--sample_name_fields` should be set to how many underscore separated parts the *sample name* has.
Note, that if the sample name had any dash (`-`) it is converted to underscore (`_`) to ensure compatibility with subsequent processes.

#### CAGEr subpipeline

To run the analysis with CAGEr with already existing bigwig files, the input is another sample sheet with 3 columns and a header row. An example is shown below (this file can be found at `docs/examples/samplesheet_cager.csv`). It is recommended to use absolute paths for the input files.

```
id,single_end,path
S10,false,[testdata/bigwig/S10.Signal.UniqueMultiple.str1.out.wig.bw testdata/bigwig/S10.Signal.UniqueMultiple.str2.out.wig.bw]
S14,false,[testdata/bigwig/S14.Signal.UniqueMultiple.str1.out.wig.bw testdata/bigwig/S14.Signal.UniqueMultiple.str2.out.wig.bw]
```

This file may be created in two steps as shown below.

```
ls /path/to/bigwig_folder/*bw > bigwigs.txt`.
./bin/sample_list_from_bigwig_list.py --filepath bigwigs.txt --singleend ["true" or "false"]
```

First, a simple txt file listing all bigwigs is created. Next, the custom script is run specifying the conversion which will create a file called `sample_list.csv`.
This script accepts two additional parameters: `--delimiter` and `--field` specifying if any part of the filename should be excluded. This might be necessary because a filename should not start with a number, as filename will be propagated as sample name and R is sensitive to names starting with a number.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run ComputationalRegulatoryGenomicsICL/customcage -params-file params.yaml -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir, default is results/)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

- [nextflow_usage](docs/nextflow_usage.md)
  - Nextflow specific information
