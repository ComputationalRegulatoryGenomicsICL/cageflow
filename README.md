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
4. Build the Bowtie2 index of the reference genome FASTA file with [`bowtie2-build`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), if the index is not provided.
5. Map the trimmed reads onto the Bowtie2 index using [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), then filter out unmapped reads and select only uniquelly mapped reads using [`samtools view`](https://www.htslib.org/doc/samtools-view.html) with options `-b -F 4 -q 20`.
6. Optionally, remove PCR and optical duplicate reads with [`samtools markdup`](https://www.htslib.org/doc/samtools-markdup.html).
7. Sort the obtained BAM files with uniquely mapped reads using [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html).
8. Index the sorted BAM files with [`samtools index`](https://www.htslib.org/doc/samtools-index.html).
9. Assess mapping quality using [`samtools stats`](https://www.htslib.org/doc/samtools-stats.html), [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html) and [`samtools idxstats`](https://www.htslib.org/doc/samtools-idxstats.html).
10. Create a CAGEexp object and call TSSs with [`CAGEr`](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) using a [BSgenome package](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) for the respective genome.
11. Create a [MultiQC](https://multiqc.info/) report.

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

The corresponding input spreadsheets can be found in [assets](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/tree/dev/assets): [samplesheet_se.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_se.csv) for single-end reads and [samplesheet_pe.csv](https://github.com/ComputationalRegulatoryGenomicsICL/customcageq/blob/dev/assets/samplesheet_pe.csv) for paired-end reads. However, you will need to use the `input_reads.sh` script to regenerate these spreadsheets with your paths to the test FASTQ files.

On these test data, CAGEr is able to call 52 TSSs with the single-end reads and 1,245 TSSs with the paired-end reads.

### How to run the pipeline

Clone the repository to your machine and use the following syntax to run the pipeline:

```bash
nextflow run customcageq/main.nf \
    --bsgenome [/path/to/]bsgenome.package[.tar.gz] \
    (--fasta /path/to/fasta/genome.fa | --index /path/to/index/bowtie2) \
    [--dedup [--dist N]] \
    --input samplesheet.csv \
    -profile <institution/docker/singularity>
```

where 
* `--bsgenome` specifies the BSgenome R package to use. If it is a file name (which should have a full path and the `.tar.gz` extension), then the package will be taken from the specified location; otherwise, the pipeline will try to install a BSgenome R package with the name `bsgenome.package` on the fly (see examples below);
* `--fasta` specifies a full path to a FASTA file containing a reference genome. This option is mandatory, unless `--index` is set. **Remark:** This option is mutually exclusive with `--index`.
* `--index` specifies a directory `bowtie2` with a Bowtie2 reference genome index. This is a mandatory option, unless `--fasta` is set. **Remark:** This option is mutually exclusive with `--fasta`.
* `--dedup` switches on PCR duplicate removal.
* `--dist N` sets an optical duplicate distance `N` to remove optical duplicates, in addition to PCR duplicates (see [`samtools markdup`](https://www.htslib.org/doc/samtools-markdup.html), option `-d`). **Remark:** The argument is optional and requires `--dedup`.
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

3. Same as example 1, but remove PCR duplicates before mapping QC and TSS calling:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer1 \
    --fasta /path/to/fasta/sacCer1.fasta \
    --dedup \
    --input customcageq/assets/samplesheet_se.csv \
    -profile docker
```

4. Same as above, but remove both PCR and optical duplicates (at a maximum distance 100, see [`samtools markdup`](https://www.htslib.org/doc/samtools-markdup.html)) before mapping QC and TSS calling:

```bash
nextflow run customcageq/main.nf \
    --bsgenome BSgenome.Scerevisiae.UCSC.sacCer1 \
    --fasta /path/to/fasta/sacCer1.fasta \
    --dedup \
    --dist 100 \
    --input customcageq/assets/samplesheet_se.csv \
    -profile docker
```

## To-do for version 2

1. Implement Damir's strategy using cutadapt to remove the first `G` and STAR for splice-aware read mapping.

2. Implement the CAGEr pipeline as a module.

3. Add plotting motifs around TSSs on both strands to check if a pyrimidine-purine (initiator-like) motif is present.

4. Check if the `nf-validation` Nextflow plugin or any other nf-core tools could help the user to create the input CSV.

5. Make it possible to run the pipeline by providing the GitHub repository name (and, possibly, a version name / commit hash), instead of making the user clone the repository first.

6. Rename `input_reads.sh` into `make_input_csv.sh` for clarity.

7. Make a "metromap" schematic of the pipeline. See, for example, the metromap for [nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.1).

8. Cite in `CITATIONS.md` all the tools that we used.

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
