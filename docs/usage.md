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
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run ComputationalRegulatoryGenomicsICL/customcage -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull ComputationalRegulatoryGenomicsICL/customcage
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [ComputationalRegulatoryGenomicsICL/customcage releases page](https://github.com/ComputationalRegulatoryGenomicsICL/customcage/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
