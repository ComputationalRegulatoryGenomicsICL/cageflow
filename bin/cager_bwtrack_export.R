save_tagclusters <- function(){
    trks <- ce |> CTSSnormalizedTpmGR("all") |> exportToTrack(ce, oneTrack = FALSE)
    for (trk in trks){
        trkname <- names(trk)
        # save tagclusters to bed
        bedfile <- file(paste0(trkname, ".bed"))
        rtracklayer::export.bed(trk, bedfile)
        # save tagclusters iq width to bed
        iqtrackfile <- file(paste0(trkname, "_iq.bed"))
        iqtrack <- exportToTrack(ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneTrack = FALSE)
        rtracklayer::export.bed(iqtrack, iqtrackfile)
        # save tagclusters to bigwig
        to_save_trk <- split(trk, strand(trk), drop = TRUE)
        bigwigfile <- file(paste0(trkname, ".bw"))
        rtracklayer::export.bw(bigwigfile)
    }
}

save_consensus_clusters <- function(){
    ctrks <- ce |> consensusClustersTpm() |> exportToTrack(ce, oneTrack = FALSE)
    # save consensus clusters to bigwig
    # save consensus clusters to bed
}
