#' Clustering tags
#'
#' @param ce normalized CAGEexp object
#' @param iqw_plot_lim Limits to IQ width plot x axis
#' @param num_core Number of cores to run on
#' @return ce clustered CAGEexp object
#' @examples
#' cager_clustering(
#' ce,
#' iqw_plot_lim = c(0, 150),
#' num_core = 4)

cager_clustering <- function(ce, iqw_plot_lim, sample_num_thr, ctss_thr, num_core){
    multicore <- TRUE
    if(num_core < 2){
        multicore <- FALSE
        num_core <- NULL
    }

    # CTSS flagging used for filtering
    ce <- CAGEr::filterLowExpCTSS(
        ce,
        thresholdIsTpm = TRUE,
        nrPassThreshold = sample_num_thr,
        threshold = ctss_thr)

    # cluster TSS with distclu
    ce <- CAGEr::distclu(
        ce,
        maxDist=20,
        keepSingletonsAbove = 5)

    # calculate IQ range
    ce <- CAGEr::cumulativeCTSSdistribution(
        ce,
        clusters = "tagClusters",
        useMulticore = multicore,
        nrCores = num_core)

    ce <- CAGEr::quantilePositions(
        ce,
        clusters = "tagClusters",
        qLow = 0.1,
        qUp = 0.9)
    
    # plot IQ width
    pdf("interquartile_width_tagclusters.pdf")
    my_plot <- plotInterquantileWidth_local(
        ce,
        clusters = "tagClusters",
        tpmThreshold = 3,
        qLow = 0.1,
        qUp = 0.9,
        xlim = iqw_plot_lim) # plot_lim
    print(my_plot)
    dev.off()

    # count CTSS
    sample_ctss_count <- list()
    for (sample in CAGEr::sampleLabels(ce)) {
        sample_ctss_count[[sample]] <- sum(as.vector(CTSStagCountDF(ce)[[sample]])>0)
    }

    sink("sample_ctss_count.txt")
    # print(data_out)
    print(sample_ctss_count)
    sink()

    plot_number_of_ctss(
        sample_ctss_count=sample_ctss_count,
        yaxistitle="Number of TSS clusters",
        mytitle="CTSS per sample",
        myfilename="ctss_counts_plots.pdf")

    # annotate tag clusters

    return(ce)
}

