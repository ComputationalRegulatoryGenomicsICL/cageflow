#' Read in Bam files to CAGEexp object
#'
#' @param bsgenome_name the name of the reference genome (bsgenome)
#' @param bam_paths list of input bam files with full path
#' @param bam_type whether it is single or Paired end
#' @param sample_names list of sample names
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
#'  bsgenome_name=hsapiens,
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

    # remove samples with empty new names
    # make it a function to call for both bams and bigwigs
    sample.idx.to.remove = which(new_names == "")
    sample_names = sample_names[-sample.idx.to.remove]
    bam_paths = bam_paths[-sample.idx.to.remove]
    bam_pairedness = bam_pairedness[-sample.idx.to.remove]

    ce = CAGEexp(genomeName     = bsgenome_name,
             inputFiles         = bam_paths,
             inputFilesType     = bam_pairedness,
             sampleLabels       = sample_names)

    ce = getCTSS(
        ce,
        removeFirstG = F,
        correctSystematicG = F,
        useMulticore = multicore,
        nrCores = cpus)

    # merge / rename samples according to user's instructioins (new_names)
    # make it a function to call for both bams and bigwigs
    name.df = data.frame(sample.name = sample_names,
                         new.name = new_names)
    
    name.df = name.df[order(name.df$new.name), ]

    name.df$merge.idx = match(name.df$new.name, unique(name.df$new.name))

    merged.sample.labels = unique(name.df$new.name)

    sample.labels = CAGEr::sampleLabels(ce)

    name.df = name.df[match(sample.labels, name.df$sample.name), ]

    merge.index = name.df$merge.idx

    ce = CAGEr::mergeSamples(ce, 
                             mergeIndex = merge.index, 
                             mergedSampleLabels = merged.sample.labels)

    return(ce)
}
