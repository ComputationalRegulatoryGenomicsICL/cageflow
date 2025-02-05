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

    # Read in TxDb object
    tx_annotation_obj <- loadDb(tx_annotation)
    ce <- annotateConsensusClusters(
        ce,
        tx_annotation_obj)
    ce <- cumulativeCTSSdistribution(
        ce,
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