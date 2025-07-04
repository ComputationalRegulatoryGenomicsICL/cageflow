#' Read in Bam files to CAGEexp object
#'
#' @param bsgenome_name the name of the reference genome (bsgenome)
#' @param bam_paths list of input bam files with full path
#' @param bam_type whether it is single or Paired end
#' @param sample_names list of sample names
#' @param new_names Character vector of new sample names (after merging/renaming).
#' @param cpus number of cores to use
#' @return a CAGEexp object
#' @examples
#' read_in_bam(
#'  bsgenome_name=hsapiens,
#'  bam_paths=["path/to/file1.bam", "path/to/file2.bam"],
#'  bam_pairedness="bam",
#'  sample_names=["S1", "S2"],
#'  cpus=1)
#' read_in_bam(
#'  bsgenome_name=BSgenome.Drerio.UCSC.danRer11,
#'  bam_paths=["path/to/file1_r1.bam", "path/to/file1_r2.bam"],
#'  bam_pairedness="bamPairedEnd",
#'  sample_names=["S1"],
#'  cpus=10)

read_in_bam <- function(
        bsgenome_name,
        bam_paths,
        bam_pairedness,
        sample_names,
        new_names,
        cpus){

    multicore <- TRUE
    if(cpus < 2){
        multicore <- FALSE
        cpus <- NULL
    }

    ce = CAGEr::CAGEexp(
        genomeName     = bsgenome_name,
        inputFiles         = bam_paths,
        inputFilesType     = bam_pairedness,
        sampleLabels       = sample_names)

    ce = CAGEr::getCTSS(
        ce,
        removeFirstG = F,
        correctSystematicG = F,
        useMulticore = multicore,
        nrCores = cpus)

    ce <- merge_labels(sample_names, new_names, ce)

    return(ce)
}
