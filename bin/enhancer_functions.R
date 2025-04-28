#' Code of exporting CAGEexp to dgCMatrix is taken from CAGEr 2.12.0
#' Enhancer calling
#'
#' @param ce initial CAGEexp object with CTSS values
#' @param cfBalanceThreshold threshold for the cagefightr balance score
#' @return enhancers
#' @examples
#' cagefightr_enhancers(
#' ce,
#' cfBalanceThreshold = 0.95,
#' )
cagefightr_enhancers <- function(
        ce,
        cfBalanceThreshold){

    # Extract CTSS count matrix as SummarizedExperiment
    se <- CAGEr::CTSStagCountSE(ce)

    # Clean up and structure metadata
    colData(se) <- colData(ce)
    rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
    colData(se)$Name <- colData(se)$sampleLabels

    # Convert counts to sparse matrix to save memory
    assays(se, withDimnames=FALSE) <- List(
        counts = as(as.matrix(as.data.frame(assays(se)[[1]])), "dgCMatrix"),
        TPM = as(as.matrix(as.data.frame(assays(se)[[2]])), "dgCMatrix"))
    
    # Save as main working object
    cfSampleCTSSs <- se

    # Calculate pooled signal across all samples (average TPM)
    cfSampleCTSSs <- CAGEfightR::calcPooled(
        cfSampleCTSSs,
        inputAssay = "TPM")
    
    # Calculate how many samples support expression at each CTSS
    cfSampleCTSSs <- CAGEfightR::calcSupport(
        cfSampleCTSSs,
        inputAssay = "counts",
        outputColumn = "support",
        unexpressed = 0)

    # Find bidirectional clusters (potential enhancers)
    sampleBCs <- CAGEfightR::clusterBidirectionally(
        cfSampleCTSSs,
        balanceThreshold = cfBalanceThreshold)

    # Filter bidirectional clusters that are supported in at least 1 sample
    finalSampleBCs <- CAGEfightR::subsetByBidirectionality(
        sampleBCs,
        samples = cfSampleCTSSs,
        minSamples = 0)

    return(finalSampleBCs)
}

#' Exclude Enhancers Overlapping Promoters
#'
#' This function filters out enhancers from CAGEfightR that overlap with consensus promoters from CAGEr
#'
#' @param BCs A Granges object representing the bidirectional clusters (BCs) from CAGEfightR.
#' @param ce A CAGEexp object from CAGEr containing the consensus clusters.
#'
#' @return A filtered version of the candidate enhancers that do not overlap with the promoters.
#' @export
#'
#' @examples
#' # Example usage:
#' # filtered_enhancers <- exclude_enhancers_overlapping_promoters(promoters, candidate_enhancers)
exclude_enhancers_overlapping_promoters <- function(BCs, ce){
    # remove parts of enhancers that overlap consensus clusters
    # true_enhancers <- GenomicRanges::setdiff(
    #     BCs,
    #     consensusClustersGR(ce),
    #     ignore.strand=TRUE)
    # OR remove regions that overlap
    promoter_overlaps <- GenomicRanges::findOverlaps(
        query=BCs,
        subject=CAGEr::consensusClustersGR(ce))
    overlapping_idx <- unique(queryHits(promoter_overlaps))
    true_enhancers <- BCs[-overlapping_idx]
    return(true_enhancers)
}

#' Annotate Enhancers with Transcript Database Information
#'
#' This function annotates a set of enhancers with information from a given transcript database.
#'
#' @param enhancers A data frame or GRanges object representing the enhancers to be annotated.
#' @param txdb A TxDb object containing transcript database information to be used for annotation.
#'
#' @return A data frame or GRanges object with the annotated enhancer information.
#' @export
#'
#' @examples
#' # Example usage:
#' # annotated_enhancers <- annotate_enhancers(enhancers, txdb)
annotate_enhancers <- function(
        enhancers,
        txdb,
        tssregion_up,
        tssregion_down){
    # TODO: maybe rather use from dELS which includes enhancers, only available for some species
    enhancerAnno <- ChIPseeker::annotatePeak(
        enhancers,
        TxDb = txdb,
        tssRegion = c(tssregion_up, tssregion_down),
        sameStrand = TRUE,
        level = "transcript",
        genomicAnnotationPriority = c(
            "Promoter", "5UTR", "3UTR",
            "Exon", "Intron",
            "Downstream", "Intergenic"))
    chipannot_plot <- ChIPseeker::plotAnnoBar(enhancerAnno)
    save_plot(
        "chipseeker_enhancer_annotation_plot.pdf",
        chipannot_plot
    )
}

#' Identify Sample-Specific Enhancers
#'
#' This function identifies sample-specific enhancers and their TPMs by comparing
#' a set of enhancers with CTSS from CAGEr
#'
#' @param true_enhancers GRanges object representing the enhancers to be compared.
#' @param ce CAGEexp object containing the sample-specific normalized CTSS data.
#'
#' @return A data frame with rows as ranges, and columns as samples, containing the TPM values for each enhancer.
#' @export
#'
identify_sample_specific_enhancers <- function(true_enhancers, ce){
    # for each sample, check the overlap between bidirectional clusters
    # and normalized CTSS counts

    # save expression per sample for each region
    enhancer_expr_per_sample <- data.frame(row.names=names(true_enhancers))

    for (sample in CAGEr::sampleLabels(ce)){
        ctss <- CAGEr::CTSSnormalizedTpmGR(ce, sample=sample)
        overlap <- GenomicRanges::findOverlaps(
            query = true_enhancers,
            subject = ctss)
        enhancer_ctss_scores <- extractList(score(ctss), overlap)
        enhancer_sample_vector <- numeric(length(enhancer_ctss_scores))
        for (i in seq_along(enhancer_ctss_scores)) {
            enhancer_sample_vector[i] <- sum(as.numeric(enhancer_ctss_scores[[i]]))
        }
        enhancer_expr_per_sample[sample] <- enhancer_sample_vector
    }
    return(enhancer_expr_per_sample)
}

#' Count Number of Enhancers Per Sample
#'
#' This function counts the number of enhancers expressed in each sample based on the enhancer expression matrix.
#'
#' @param enhancer_expr_per_sample A data frame where rows represent enhancers and columns represent samples, containing expression values.
#'
#' @return A named list where each element corresponds to a sample and contains the count of expressed enhancers.
#' @export
#'
#' @examples
#' # Example usage:
#' # enhancer_counts <- count_number_of_enhancers(enhancer_expr_per_sample)
count_number_of_enhancers <- function(enhancer_expr_per_sample) {
    sample_enhancer_count <- list()
    for (sample in colnames(enhancer_expr_per_sample)) {
        sample_enhancer_count[[sample]] <- sum(as.vector(enhancer_expr_per_sample[,sample]) > 0)
    }
    return(sample_enhancer_count)
}

