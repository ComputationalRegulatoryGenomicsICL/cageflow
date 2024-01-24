# &beta;

## Introduction

**ComputationalRegulatoryGenomicsICL/customcageq** is a Nextflow pipeline to process CAGE sequencing data from raw reads to the creation of a CAGEexp (CAGEr) object containing called TSSs. The pipeline is specifically designed to be used upstream of CAGEr.

### Input

Either single-end or paired-end raw CAGE reads. Only one type of reads (either single- or paired-end) can be used in one run of the pipeline.

### Output

A CAGEexp (CAGEr) object with called TSSs, ready for a downstream analysis with CAGEr. The intermediate and final results are stored in the `results` directory. The final CAGEexp object is stored in an RDS file in the `results/cager` directory.

### Steps

1. Merge per-lane FASTQ files with the [`nf-core/cat_fastq`](https://nf-co.re/modules/cat_fastq) module.
2. Report raw read quality with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
3. Trim adapters with [`TrimGalore`](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) and run [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on trimmed reads.
4. Download the reference genome FASTA file from the UCSC server, if not provided locally, using [`BuxyBox wget`](https://boxmatrix.info/wiki/Property:wget) within a custom module [`DOWNLOAD_FASTA`](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/modules/local/downloadfasta.nf).
5. Build the Bowtie2 index of the reference genome FASTA file with [`bowtie2-build`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), if the index is not provided locally.
6. Map the trimmed reads onto the Bowtie2 index using [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), then filter out unmapped reads and select only uniquelly mapped reads using [`samtools view`](https://www.htslib.org/doc/samtools.html) with options `-b -F 4 -q 20`.
7. Sort the obtained BAM files with uniquelly mapped reads using [`samtools sort`](https://www.htslib.org/doc/samtools.html).
8. Index the sorted BAM files with [`samtools index`](https://www.htslib.org/doc/samtools.html).
9. Create a CAGEexp object and call TSSs with [`CAGEr`](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) using a [BSgenome package](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) for the respective genome.

## Usage

### Prepare for your first run

Currently, pipeline works with Nextflow v23.04. Make sure that you have the latest version of Docker (if running the pipeline on a laptop / PC) or Singularity (if running on a high-performance cluster).

### Prepare your input data

Prepare the sample sheet with the description of input samples. In case of single-end reads, it should look like this:

```csv
sample,fastq_1,fastq_2,single_end
S1,/path/to/fastq/S1_S1_L001_R1_001.fastq.gz,,True
S1,/path/to/fastq/S1_S1_L002_R1_001.fastq.gz,,True
S2,/path/to/fastq/S2_S2_L001_R1_001.fastq.gz,,True
S2,/path/to/fastq/S2_S2_L002_R1_001.fastq.gz,,True
```

where
* `sample` is a unique identifier of a sample;
* `fastq_1` (and `fastq_2` in the case of paired-end reads) is a full path to the read libraries. In case of paired-end reads, `fastq_1` contains the full path to forward reads, while `fastq_2` contains the full path to reverse reads. One sample can be represented by more than one library if lanes are stored separately;
* `single_end` should be set to `True` for single-end reads and to `False` for paired-end reads.

For paired-end reads, `fastq_2` should contain the full path to reverse reads, while `single_end` should be set to `False`.

You can generate the input CSV table automatically using the [`input_reads.sh`](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/bin/input_reads.sh) script. It takes two positional arguments:

```bash
input_reads.sh /path/to/fastq_dir /path/to/samplesheet.csv
```

where
* `/path/to/fastq_dir` is a full path to a directory with raw FASTQ files;
* `/path/to/samplesheet.csv` is a file name, with a full path, of a CSV file to create.

To run the script as a standalone executable (that is, without the need to write `bash` before its name), add execution permissions to the script after cloning the repository:

```bash
chmod +x input_reads.sh
```

### Toy input data for testing

The pipeline has toy *S. cerevisiae* CAGE data stored in [assets/sacer_fq](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq) for testing purposes (single-end reads in the [se](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq/se) subfolder and paired-end reads in [pe](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq/pe) subfolder). The single-end reads were obtained by subsampling the [ERR2495152](https://www.ebi.ac.uk/ena/browser/view/ERR2495152) dataset published by ([Börlin et al., 2018](https://academic.oup.com/femsyr/article/19/2/foy128/5257840)), while the paired-end reads were obtained by subsampling the [SRR1631657](https://www.ebi.ac.uk/ena/browser/view/SRR1631657) dataset published by ([Chabbert et al., 2015](https://www.embopress.org/doi/full/10.15252/msb.20145776)).

The corresponding input spreadsheets can be found in [assets](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets): [samplesheet_se.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_se.csv) for single-end reads and [samplesheet_pe.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_pe.csv) for paired-end reads.

On these test data, CAGEr is able to call 52 TSSs with the single-end reads and 1,245 TSSs with the paired-end reads.

### How to run the pipeline

Clone the repository to your machine and use the following syntax to run the pipeline:

```bash
nextflow run customcageq/main.nf \
    --bsgenome [/path/to/]bsgenome.package[.tar.gz] \
    [--fasta /path/to/fasta/genome.fa | --index /path/to/index/genome] \
    --input samplesheet.csv \
    -profile <institution/docker/singularity>
```

where 
* `--bsgenome` specifies the BSgenome R package to use. If it is a file name (which should have a full path and the `.tar.gz` extension), then the package will be taken from the specified location; otherwise, the pipeline will try to install a BSgenome R package with the name `bsgenome.package` on the fly (see examples below);
* `--fasta` specifies a full path to a FASTA file containing the reference genome. This is an optional parameter. If it is specified, then the pipeline will take the reference genome FASTA from `/path/to/fasta/genome.fa` and use it to create the Bowtie2 index. **Remark:** This option is mutually exclusive with `--index`.
* `--index` specifies a directory with a Bowtie2 reference genome index. This is an optional parameter. If it is specified, then the pipeline will use this index to map the input CAGE data. **Remark:** This option is mutually exclusive with `--fasta`.
* `--input` specifies the input CSV samplesheet.
* `-profile` is a Nextflow option that specifies a config file to use with Nextflow on a given machine. See [`nf-core/configs`](https://github.com/nf-core/configs) for ready-to-use institutional configs, including the one for Jex (the high-performance computing cluster of the [Laboratory of Medical Sciences](https://lms.mrc.ac.uk/)). Also, see the [Jex wiki](https://hpcwiki.lms.mrc.ac.uk/docs/software/software/workflow_managers/#nextflow) on how to run Nextflow on Jex. Alternatively, this option can be used to specify the containerization technology to use.

### Examples

1. Call TSSs from the test yeast single-end CAGE reads using a locally stored reference FASTA file and the `BSgenome.Scerevisiae.UCSC.sacCer1` R package. The package is automatically installed within the CAGEr container on the fly and used there with CAGEr:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer1 \
    --fasta /path/to/fasta/sacCer1.fasta \
    --input customcageq/assets/samplesheet_se.csv \
    -profile docker
```

2. Call TSSs from the test yeast paired-end CAGE reads using a locally stored Bowtie2 index and the locally stored `BSgenome.Scerevisiae.UCSC.sacCer1` R package. The package is automatically installed from the `.tar.gz` archive within the CAGEr container and used with CAGEr:

```bash
nextflow run customcageq/main.nf \
    --bsgenome /path/to/bsgenome/BSgenome.Scerevisiae.UCSC.sacCer1_1.4.0.tar.gz \
    --index /path/to/index/bowtie2 \
    --input customcageq/assets/samplesheet_pe.csv \
    -profile docker
```

3. Call TSSs from the test yeast single-end CAGE reads using the [sacCer1 reference FASTA file](http://hgdownload.cse.ucsc.edu/goldenPath/sacCer1/bigZips/sacCer1.fa.gz) automatically downloaded from the UCSC server and the `BSgenome.Scerevisiae.UCSC.sacCer1` R package automatically installed within the CAGEr container on the fly and used there with CAGEr. As neither `--fasta`, nor `--index`, are specified, the genome name `sacCer1` is automatically taken from the provided BSgenome package name (or the file name of the package `.tar.gz` archive stored locally, if provided) to find and download the respective reference FASTA file:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer1 \
    --input customcageq/assets/samplesheet_se.csv \
    -profile docker
```

## Fix in &beta;

1. **[Done]** Support both Docker and Singulatiry for the CAGEr container. The problem is that the module cannot install a required BSgenome if run in Singularity.

2. **[Done]** Describe in this README how to use the `input_reads.sh` script.

3. **[Done]** Adjust default resource allocation for a generic HPC.

## To-do for version 2

1. [**In progress**] Make the pipeline compatible with Nextflow v23.10.0 (or later). The problem with this version of Nextflow is that the reference genome index is not replicated in the corresponding input channel of the nf-core module bowtie2align according to the number of samples to map. Therefore, only one sample gets mapped.

2. Add mapping stats to the MultiQC report. For this, use nf-core modules [samtools/flagstat](https://nf-co.re/modules/samtools_flagstat) to count the number of alignments for each FLAG type, [samtools/idxstats](https://nf-co.re/modules/samtools_idxstats) to print mapping stats per chromosome and [samtools/stats](https://nf-co.re/modules/samtools_stats) to print a comprehensive mapping summary.

3. Add plotting motifs around TSSs on both strands to check if a pyrimidine-purine (initiator-like) motif is present.

4. Check if the `nf-validation` Nextflow plugin or any other nf-core tools could help the user to create the input CSV.

5. Make it possible to run the pipeline by providing the GitHub repository name (and, possibly, a version name / commit hash), instead of making the user clone the repository first.

6. Move the code base into a current nf-core template (so that anyone forking our repo could use nf-core tools to develop their derivative pipeline).

7. Rename `input_reads.sh` into `make_input_csv.sh` for clarity.

8. Make a "metromap" schematic of the pipeline. See, for example, the metromap for [nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.1).

9. Cite in `CITATIONS.md` all the tools that we used.

## Possible future features

1. Improve the `input_reads.sh` script according to Damir's comments that he left within it. 

2. Use the iGenomes repository for obtaining reference genomes: we could obtain ready-made Bowtie2 indexes from there, instead of FASTA files (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/README.txt). But we need to make sure that for a given reference genome we download only the Bowtie2 index and not the whole bundle of files (which also includes other indexes and a genome sequence). iGenomes, in general, may be tricky to use because the genomes there may be outdated or have other problems (see the notifications in the nf-core documentation: https://nf-co.re/docs/usage/reference_genomes), so I would not prioritise the use of iGenomes, although being able to download a ready-made Bowtie2 index would be very useful!

3. Implement BSgenome forging based on a seed file specified by the user. Apart from other fields, a seed file contains a path to the directory with a FASTA file or a 2bit file to forge the BSgenome (source directory). Provide for amending the path to the source directory using an optional `--sourcedir` parameter. Additionally, provide for the BSgenome forging and CAGE data preprocessing in one go with an option `--forge` and for BSgenome forging only - with an option `--forge-only`. 

* Forge a BSgenome and proceed with CAGE preprocessing using the forged BSgenome and a FASTA file. The FASTA file should either be provided for the BSgenome forging and, hence, located in the `source_dir`, or be provided with the `--fasta` option:

```bash
--forge \
--bsgenome mySpecies \
[--sourcedir /path/to/source_dir] \
[--fasta /path/to/fasta/genome.fasta] \
--seed /path/to/seed/mySpecies_seed.txt \
--input input.csv
```

A FASTA file must be specified with the `--fasta` option if it was not provided for the BSgenome forging.

* Only forge a BSgenome and exit:

```bash
--forge-only \
--bsgenome mySpecies \
[--sourcedir /path/to/source_dir] \
--seed /path/to/seed/mySpecies_seed.txt
```

The forged package would be called `BSgenome.LatinName.custom.mySpecies`, where `LatinName` could be something like `Scerevisia` (taken from the seed file) and `mySpecies` could be something like `sacCer2`.

## Credits

**ComputationalRegulatoryGenomicsICL/customcageq** was originally written by Pavel Nikitin ([@nikitin-p](https://github.com/nikitin-p)), Sviatoslav Sidorov ([@sidorov-si](https://github.com/sidorov-si)) and Damir Baranasic ([@da-bar](https://github.com/da-bar)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  ComputationalRegulatoryGenomicsICL/customcage for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
