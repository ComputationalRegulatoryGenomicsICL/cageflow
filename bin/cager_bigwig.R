# Function from CAGEr version 2.12.0
#' @name getRefGenome
#' 
#' @title Return a BSgenome or throws an error
#' 
#' @details If the `reference.genome` object exists and is a BSgenome_, it will
#' be returned.  This allows the user to run things like
#' `seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg19) <- "NCBI"` on the BSgenome
#' object before running `getCTSS`.  If the `reference.genome` object does not
#' exist, attempts to load it and return it, or throws an error if not available.
#' 
#' @return A BSgenome object of the same name as the `reference.genome` argument.
#' 
#' @param reference.genome
#' 
#' @author Charles Plessy
#' @importFrom utils installed.packages
#' 
#' @noRd

getRefGenome <- function(reference.genome) {
  if (is.null(reference.genome))
    stop("Can not run this function with a NULL genome; see ", sQuote('help("genomeName")'), ".")
  if (exists(reference.genome))
    if ("BSgenome" %in% class(get(reference.genome)))
      return(get(reference.genome))
    else
      stop("The ", sQuote(reference.genome), " object in the namespace is not a BSgenome")
  if(reference.genome %in% rownames(installed.packages()) == FALSE)
    stop("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
  requireNamespace(reference.genome)
  getExportedValue(reference.genome, reference.genome)
}

# Function from CAGEr version 2.12.0
#' coerceInBSgenome
#' 
#' A private (non-exported) function to discard any range that is
#' not compatible with the CAGEr object's BSgenome.
#' 
#' @param gr The genomic ranges to coerce.
#' @param genome The name of a BSgenome package, which must me installed,
#'   or \code{NULL} to skip coercion.
#' 
#' @return A GRanges object in which every range is guaranteed to be compatible
#' with the given BSgenome object.  The sequnames of the GRanges are also set
#' accordingly to the BSgenome.
#' 
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors %in%

coerceInBSgenome <- function(gr, genome) {
  if (is.null(genome)) return(gr)
  genome <- getRefGenome(genome)
  gr <- gr[seqnames(gr) %in% seqnames(genome)]
  gr <- gr[! end(gr) > seqlengths(genome)[as.character(seqnames(gr))]]
  seqlevels(gr) <- seqlevels(genome)
  seqinfo(gr) <- seqinfo(genome)
  gr
}

#' Read in BigWig files to CAGEexp object
#'
#' @param bsgenome_name the name of the reference genome (bsgenome)
#' @param bigwig_paths list of input bigwig files with full path
#' @param sample_names_files_dict dictionary of matching sample names to files
#' @return a CAGEexp object
#' @examples
#' read_in_bigwig(
#'  bsgenome_name="BSgenome.Scerevisiae.UCSC.sacCer3",
#'  bigwig_paths=["path/to/file1.bw", "path/to/file2.bw"]
#'  )

read_in_bigwig <- function(
    bsgenome_name,
    bigwig_paths,
    sample_names_files_dict){
    # action){

  bigwig_paths <- stringr::str_squish(bigwig_paths)

  bigwigs = unlist(
    stringr::str_split(
      stringr::str_remove_all(
        bigwig_paths, ","),
      stringr::fixed(" ")))

  signals = lapply(
    bigwigs,
    function(x) {
      track_in <- rtracklayer::import(x)
      track_bs <- coerceInBSgenome(track_in, bsgenome_name)
    })

  signal_names <- c()
  for (bn in basename(bigwigs)){
    signal_names <- append(signal_names, sample_names_files_dict[[bn]])
  }
  names(signals) = signal_names

  # print(signals)
  # print(grepl("str1", names(signals)))

  signalsSplit = split(
    signals,
    grepl("str1", names(signals)))

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
  } else if (grepl( "_str1", names(plus)[1], fixed = TRUE)) {
    plus_sample_names = stringr::str_remove_all(
      names(plus),
      "_str1")
    minus_sample_names = stringr::str_remove_all(
      names(minus),
      "_str2")
  } else {
    plus_sample_names = names(plus)
    minus_sample_names = names(minus)
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

  # Step 0: Create a CAGEexp object, filenames only of str1
  ce <- new( "CAGEexp"
     , colData = DataFrame( inputFiles     = bigwigs[grep("str1", bigwigs)]
                          , sampleLabels   = plus_sample_names
                          , inputFilesType = "CTSStable"
                          , row.names      = plus_sample_names)
     , metadata = list(genomeName = bsgenome_name))

  # Step 1: Load each file as GRangesList where each GRange is a CTSS data.
  merged = mapply(c, plus, minus)
  merged_gpos <- lapply(merged, function(x) {
    gp <- GPos(stitch=FALSE, x)
    score(gp) <- x$score
    gp <- coerceInBSgenome(gp, bsgenome_name)
    gp <- sort(gp)
  })
  merged_gpos <-  GRangesList(merged_gpos)

  # Step 2: Create GPos representing all the nucleotides with CAGE counts in the list.
  rowRanges <- sort(unique(unlist(merged_gpos)))
  mcols(rowRanges) <- NULL

  # Step 3: Fold the GRangesList in a expression DataFrame of Rle-encoded counts.
  assay <- DataFrame(V1 = Rle(rep(0L, length(rowRanges))))
  expandRange <- function(global, local) {
    x <- Rle(rep(0L, length(global)))
    x[global %in% local] <- score(local)
    x
  }
  for (i in seq_along(merged_gpos))
    assay[,i] <- expandRange(rowRanges, merged_gpos[[i]])
  
  rowRanges <- new("CTSS", rowRanges, bsgenomeName = bsgenome_name)
  colnames(assay) <- names(merged)

  # Setp 4: Put the data in the appropriate slot of the MultiAssayExperiment.
  CTSStagCountSE(ce) <- SummarizedExperiment(
      rowRanges = rowRanges,
      assays    = SimpleList(counts = assay))
  
  # Step 5: update the sample metadata (colData).
  ce$librarySizes <- unlist(lapply(CTSStagCountDF(ce), sum))
  
  # Merge if necessary
  # if (action is not NULL){
  #     ce <- mergeSamples(ce, mergeIndex = seq(1,len(camplename)),#c(3,2,4,4,1), 
  #                 mergedSampleLabels = action)#c("Zf.unfertilized.egg", "Zf.high", "Zf.30p.dome", "Zf.prim6"))
  # }

  # Setp 6: Return the modified object.
  ce
}
