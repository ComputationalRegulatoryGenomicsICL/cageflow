#' Read in Bam files to CAGEexp object
#'
#' @param reference_name the name of the reference genome (bsgenome)
#' @param input_files list of input bam files with full path
#' @param bam_type whether it is single or Paired end
#' @param sample_names list of sample names
#' @param cpus number of cores to use
#' @return a CAGEexp object
#' @examples
#' read_in_bam(
#'  bsgenome_name=hsapiens,
#'  input_files=[path/to/file1, path/to/file2],
#'  bam_pairedness=bam,
#'  sample_names=[S1, S2],
#'  cpus=1)
#' read_in_bam(
#'  bsgenome_name=hsapiens,
#'  input_files=[path/to/file1_r1, path/to/file1_r2],
#'  bam_pairedness=bamPairedEnd,
#'  sample_names=[S1],
#'  cpus=10)

read_in_bam <- function(
        bsgenome_name,
        input_files,
        bam_pairedness,
        sample_names,
        cpus){

    ce = CAGEexp(genomeName     = bsgenome_name,
             inputFiles     = input_files,
             inputFilesType = bam_pairedness,
             sampleLabels   = sample_names)

    ce = getCTSS(
        ce,
        removeFirstG = F,
        correctSystematicG = F,
        useMulticore = T,
        nrCores = cpus)

    return(ce)
}

