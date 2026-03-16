#' Extract CTSS score per sample per region

overlapper <- function(ce, sample, remove_gg_initiator, keep_only_yr_yc){
    ctss_raw <- CAGEr::CTSSnormalizedTpmGR(ce, sample=sample)
    # TODO: this code is kind of repeated 3 times: for tagclustering, enhancer calling and consensus cluster scoring
    if (remove_gg_initiator) {
        dinuc <- ctss_raw %>%
            GRanges() %>%
            promoters(upstream = 1, downstream = 1) %>%
            {getSeq(BSgenome.Hsapiens.UCSC.hg38, trim(.))}
            # {getSeq(BSgenome.Drerio.UCSC.danRer11, trim(.))}
        ctss_raw$dinuc <- as.character(dinuc)
        # when GG is the starting dinucleotide, the flag is set to FALSE
        not_gg_start <- !(ctss_raw$dinuc == "GG")
        # do not filter when dinucleotide is missing
        not_gg_start[is.na(not_gg_start)] <- TRUE
        # Apply GG initiation filter
        ctss <- ctss_raw[not_gg_start]
    } else if (keep_only_yr_yc) {
        dinuc <- ctss_raw %>%
            GRanges() %>%
            promoters(upstream = 1, downstream = 1) %>%
            {getSeq(BSgenome.Hsapiens.UCSC.hg38, trim(.))}
            # {getSeq(BSgenome.Drerio.UCSC.danRer11, trim(.))}
        ctss_raw$dinuc <- as.character(dinuc)
        # in YR group, select: CG, CA, TG, TA
        # in YC group, select: CC, TC
        yryc_group <- c("CG", "CA", "TG", "TA", "CC", "TC")
        # when YR or YC is the starting dinucleotide, the flag is set to TRUE, else FALSE
        yr_yc_start <- (ctss_raw$dinuc %in% yryc_group)
        # Apply initiation filter
        ctss <- ctss_raw[yr_yc_start]
    } else {
        ctss <- ctss_raw
    }

    overlap <- GenomicRanges::findOverlaps(
            query = consensusClustersGR(ce),
            subject = ctss)
    grouped_scores <- extractList(score(ctss), overlap)
    return(grouped_scores)
}

sum_scores_per_location <- function(scores){
    # sum values
    result_vector <- numeric(length(scores))
    # Loop through each element and add the value at the correct position
    for (i in seq_along(scores)) {
        result_vector[i] <- sum(as.numeric(scores[[i]]))
    }
    return(result_vector)
}

# Replace scores in consensus clusters with the sum of normalized CTSS scores
score_from_ctss <- function(ce, remove_gg_initiator, keep_only_yr_yc){
    # use normalized CTSS to calculate and update score
    ctss_score_df <- data.frame(names = names(consensusClustersGR(ce)))
    sample_sum <- 0

    for (sample in CAGEr::sampleLabels(ce)){
        new_scores <- overlapper(ce, sample, remove_gg_initiator, keep_only_yr_yc)
        ctss_score_df[[sample]] <- sum_scores_per_location(new_scores)
        sample_sum <- sample_sum + ctss_score_df[[sample]]
    }

    write.csv(
        ctss_score_df,
        "tables/consensus_cluster_per_sample_ctss.csv",
        row.names = FALSE,
        quote = FALSE)

    # Convert to Rle and assign to consensus cluster
    # add to a new column called ctss_score
    elementMetadata(consensusClustersGR(ce))[["ctss_score"]] <- Rle(sample_sum)

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
        iqhigh,
        remove_gg_initiator,
        keep_only_yr_yc){

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

    ce <- score_from_ctss(ce, remove_gg_initiator, keep_only_yr_yc)

    # Read in TxDb object
    tx_annotation_obj <- loadDb(tx_annotation)
    ce <- annotateConsensusClusters(
        ce,
        tx_annotation_obj)

    ce <- CAGEr::cumulativeCTSSdistribution(
        ce,
        clusters = "consensusClusters",
        useMulticore = multicore,
        nrCores = num_core)
    ce <- CAGEr::quantilePositions(
        ce,
        clusters = "consensusClusters",
        qLow = iqlow,
        qUp = iqhigh,
        useMulticore = multicore,
        nrCores = num_core)

    return(ce)
}