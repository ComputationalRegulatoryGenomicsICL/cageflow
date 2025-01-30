bigwig_export <- function(x, y){
    tracks <- split(x, strand(x))
    rtracklayer::export.bw(tracks$`+`, paste0(y,"_plus.bw"))
    min <- tracks$`-`
    min$score <- min$score * (-1)
    rtracklayer::export.bw(tracks$`+`, paste0(y,"_plus.bw"))
    rtracklayer::export.bw(min, paste0(y, "_minus.bw") )
}

export_tagclusters <- function(ce){
    # export normalized TSS counts into bigwig
    mapply(
        bigwig_export,
        CAGEr::CTSSnormalizedTpmGR(ce, "all"),
        CAGEr::sampleLabels(ce))

    # export raw TSS counts into bigwig
    mapply(
        bigwig_export,
        CAGEr::CTSStagCountGR(ce, "all"),
        CAGEr::sampleLabels(ce))

    bedTracks <- CAGEr::exportToTrack(
        ce, 
        what = "tagClusters", 
        qLow = 0.1, qUp = 0.9, 
        oneTrack = FALSE)

    mapply(function(x, y){
        rtracklayer::export.bed(x, paste0(y,"_tagClusters.bed"))
    }, bedTracks, CAGEr::sampleLabels(ce))
}

export_consensus_clusters <- function(ce){
    ccbedTracks <- CAGEr::exportToTrack(
        ce, 
        what = "consensusClusters", 
        colorByExpressionProfile = TRUE,
        oneTrack = FALSE)

    rtracklayer::export.bed(ccbedTracks, "consensusClusters.bed")
}
