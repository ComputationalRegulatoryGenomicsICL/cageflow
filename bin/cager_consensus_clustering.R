#' Extract CTSS score per sample per region

overlapper <- function(ce, sample){
    ctss <- CAGEr::CTSSnormalizedTpmGR(ce, sample=sample)
    overlap <- GenomicRanges::findOverlaps(
            query = consensusClustersGR(ce),
            subject = ctss)
    grouped_scores <- extractList(score(ctss), overlap)
    return(grouped_scores)
}

# Replace scores in consensus clusters with the sum of normalized CTSS scores
score_from_ctss <- function(ce){
    # use normalized CTSS to calculate and update score
    # extract the scores of the first sample
    sample <- CAGEr::sampleLabels(ce)[[1]]
    grouped_scores <- overlapper(ce, sample)

    # ...and of the rest of the samples
    for (sample in CAGEr::sampleLabels(ce)[2:length(CAGEr::sampleLabels(ce))]){
        new_scores <- overlapper(ce, sample)
        combined <- mapply(append, grouped_scores, new_scores, SIMPLIFY=FALSE)
        grouped_scores <- RleList(combined)
    }

    # sum values
    result_vector <- numeric(length(grouped_scores))
    # Loop through each element and add the value at the correct position
    for (i in seq_along(grouped_scores)) {
        result_vector[i] <- sum(as.numeric(grouped_scores[[i]]))
    }

    # Convert to Rle  and assign to consensus cluster
    # TODO: this does not actually assign the value which is 
    # interesting because within cager it does
    score(consensusClustersGR(ce)) <- Rle(result_vector)

    return(ce)
}

#' Consensus clustering of CTSS across samples
#'
#' @param ce normalized CAGEexp object
#' @param tpmThreshold Threshold for calling consensus clusters (in aggregateTagClusters function)
#' @param maxDist Maximum distance for calling consensus clusters (in aggregateTagClusters function)
#' @param tx_annotation SQLite file with a TxDb genome annotation package
#' @return ce clustered CAGEexp object
#' @examples
#' consensus_clustering(
#' ce,
#' tpmThreshold = 5,
#' maxDist = 100,
#' tx_annotation = "annotation_from_gtf.sqlite",
#' num_core = 4)

consensus_clustering <- function(
        ce,
        tpmThreshold,
        maxDist,
        tx_annotation,
        num_core,
        iqlow,
        iqhigh){

    multicore <- TRUE
    if(num_core < 2){
        multicore <- FALSE
        num_core <- NULL
    }

    ce <- CAGEr::aggregateTagClusters(
        ce,
        tpmThreshold = tpmThreshold,
        qLow = iqlow,
        qUp = iqhigh,
        maxDist = maxDist)

    ce <- score_from_ctss(ce)

    # Read in TxDb object
    tx_annotation_obj <- loadDb(tx_annotation)
    ce <- annotateConsensusClusters(
        ce,
        tx_annotation_obj)
    ce <- cumulativeCTSSdistribution(
    clusters = "consensusClusters",
        useMulticore = multicore,
        nrCores = num_core)
    ce <- quantilePositions(
        ce,
        clusters = "consensusClusters",
        qLow = iqlow,
        qUp = iqhigh,
        useMulticore = multicore,
        nrCores = num_core)

    return(ce)
}