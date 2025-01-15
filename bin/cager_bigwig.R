#' Read in BigWig files to CAGEexp object
#'
#' @param bsgenome_name the name of the reference genome (bsgenome)
#' @param bigwig_paths list of input bigwig files with full path
#' @param chromosome_names list of chromosome names to keep
#' @param cpus number of cores to use
#' @return a CAGEexp object
#' @examples
#' read_in_bigwig(
#'  bsgenome_name="BSgenome.Scerevisiae.UCSC.sacCer3",
#'  bigwig_paths=["path/to/file1.bw", "path/to/file2.bw"],
#'  cpus=4)

read_in_bigwig <- function(
    bsgenome_name,
    bigwig_paths,
    chromosome_names,
    cpus){

  bigwig_paths <- stringr::str_squish(bigwig_paths)

  bigwigs = unlist(
    stringr::str_split(
      stringr::str_remove_all(
        bigwig_paths, ","),
      stringr::fixed(" ")))

  signals = lapply(
    bigwigs,
    function(x) rtracklayer::import(x))

  names(signals) = basename(bigwigs)

  signals_chr_filt <- list()
  for (sname in names(signals)){
    signal <- signals[[sname]]
    signal <- signal[seqnames(signal) %in% chromosome_names]
    signals_chr_filt[[sname]] <- signal
  }
  

  signalsSplit = split(
    signals_chr_filt,
    grepl("str1", names(signals_chr_filt)))

  plus = lapply(signalsSplit$`TRUE`, function(x) {
    strand(x) = "+"
    return(x)
  })

  minus = lapply(signalsSplit$`FALSE`, function(x) {
    strand(x) = "-"
    return(x)
  })

  if (grepl( ".Signal.UniqueMultiple.str1.out.wig.bw", names(plus)[1], fixed = TRUE)){
    plus_sample_names = stringr::str_remove_all(
      names(plus),
      ".Signal.UniqueMultiple.str1.out.wig.bw")
    minus_sample_names = stringr::str_remove_all(
      names(minus),
      ".Signal.UniqueMultiple.str2.out.wig.bw" )
  } else if (grepl( ".Signal.Unique.str1.out.wig.bw", names(plus)[1], fixed = TRUE)) {
    plus_sample_names = stringr::str_remove_all(
      names(plus),
      ".Signal.Unique.str1.out.wig.bw")
    minus_sample_names = stringr::str_remove_all(
      names(minus),
      ".Signal.Unique.str2.out.wig.bw")
  }

  names(plus) <- plus_sample_names
  names(minus) <- minus_sample_names

  if (!all(plus_sample_names == minus_sample_names)) {
    if (setequal(plus_sample_names, minus_sample_names)) {
      minus = minus[match(plus_sample_names, minus_sample_names)]
    } else {
      stop("Error: Some basenames of minus- and plus-strand bigWigs are different! Are these bigWigs from different sets of samples? Exit.")
    }
  }

  merged = mapply(c, plus, minus)

  ctss.obj = purrr::reduce(
    merged,
    function(x, y) {
      dplyr::full_join(
        as.data.frame(x),
        as.data.frame(y),
        by = c("seqnames", "start", "strand")) %>%
      dplyr::select(
        c("seqnames", "start", "strand"),
        contains("score"))
    }) %>%
    dplyr::rename(
      chr = "seqnames",
      pos = "start") %>%
    dplyr::mutate(
      across(
        where(is.numeric),
        ~tidyr::replace_na(., 0)))

  names(ctss.obj) = c("chr", "pos", "strand", names(merged))

  ctss.obj %<>%
    dplyr::mutate(across(all_of(plus_sample_names), as.integer))
  ctss.obj$chr = as.character(ctss.obj$chr)
  ctss.obj$strand = as.character(ctss.obj$strand)
  cageexpobj = as(ctss.obj, "CAGEexp")

  colData(cageexpobj)$inputFilesType <- "CTSStable"

  metadata(cageexpobj)$genomeName = bsgenome_name

  rowRanges(cageexpobj@ExperimentList$tagCountMatrix) = as(
    rowRanges(cageexpobj@ExperimentList$tagCountMatrix),
    Class = "CTSS")
  
  return(cageexpobj)
}
