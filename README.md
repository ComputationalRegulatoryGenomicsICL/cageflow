## Introduction

**ComputationalRegulatoryGenomicsICL/customcageq** is a Nextflow pipeline to process CAGE sequencing data from raw reads to the creation of a CAGEexp (CAGEr) object containing called TSSs. The pipeline is specifically designed to be used upstream of CAGEr.

### Input

Either single-end or paired-end raw CAGE reads. Onle one type of reads (either single- or paired-end) can be used in one run of the pipeline.

### Output

A CAGEexp (CAGEr) object with called TSSs, ready for a downstream analysis with CAGEr.

### Steps

1. Merge per-lane FASTQ files with the [`nf-core/cat_fastq`](https://nf-co.re/modules/cat_fastq) module.
2. Report raw read quality with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
3. Trim adapters with [`TrimGalore`](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) and run [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on trimmed reads.
4. Download the reference genome FASTA file from UCSC, if not provided locally, using [`BuxyBox wget`](https://boxmatrix.info/wiki/Property:wget) within a custom module [`DOWNLOAD_FASTA`](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/modules/local/downloadfasta.nf).
5. Build the Bowtie2 index of the reference genome FASTA file with [`bowtie2-build`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), if the index is not provided locally.
6. Map the trimmed reads onto the Bowtie2 index using [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) with options `-b -F 4 -q 20` to filter out unmapped reads and select only uniquelly mapped reads.
7. Sort the obtained BAM files with uniquelly mapped reads using [`samtools sort`](https://www.htslib.org/doc/samtools.html).
8. Index the sorted BAM files with [`samtools index`](https://www.htslib.org/doc/samtools.html).
9. Create a CAGEexp object and call TSSs with [`CAGEr`](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) using a [BSgenome package](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) for the respective genome.

## Usage

### Prepare for your first run

Make sure you have the latest version of Nextflow, as well as the latest version of Docker (if running the pipeline on a laptop / PC) or Singularity (if running on a high-performance cluster).

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
* `fastq_1` (and `fastq_2` in the case of paired-end reads) is a full path to the read libraries. In case of paired-end reads, `fastq_1` contains the full path to forward reads, while `fastq_2` contains the full path to reverse reads. One sample can be represented by more than one library if each lane stored separately;
* `single_end` should be set to `True` for single-end reads and to `False` for paired-end reads.

For paired-end reads, `fastq_2` should contain the full path to reverse reads, while `single_end` should be set to `False`.

You can generate the input CSV table automatically using the [`input_reads.sh`](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/bin/input_reads.sh) script.

### Toy input data for testing

The pipeline has toy *S. cerevisiae* "CAGE" data stored in [assets/sacer_fq](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq) for testing purposes (single-end reads in the [se](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq/se) subfolder and paired-end reads in [pe](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets/sacer_fq/pe) subfolder). The data was obtained from the *S. cerevisiae* genome bioinformatically, by random sampling of its subsequences.

The corresponding input spreadsheets can be found in [assets](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets): [samplesheet_se.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_se.csv) for single-end reads, and [samplesheet_pe.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_pe.csv) for paired-end reads.

On these data, CAGEr is able to call several tens of TSSs.

### Synopsis

To run the pipeline, use the following syntax

```bash
nextflow run ComputationalRegulatoryGenomicsICL/customcage \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv
```

### Run examples

...

## Credits

ComputationalRegulatoryGenomicsICL/customcage was originally written by Damir Baranasic.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  ComputationalRegulatoryGenomicsICL/customcage for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
