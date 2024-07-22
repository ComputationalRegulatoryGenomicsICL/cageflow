#' Clustering tags
#'
#' @param ce normalized CAGEexp object
#' @param threshold
#' @param thresholdIsTpm
#' @param nrPassThreshold
#' @param method
#' @param maxDist
#' @param removeSingletons
#' @param keepSingletonsAbove
#' @param qLow
#' @param qUp
#' @param tpmThreshold
#' @param plot_lim
#' @param num_core
#' @param cagerFolder
#' @return ce clustered CAGEexp object
#' @examples
#' cager_clustering(
#' threshold = 1,
#' thresholdIsTpm = TRUE,
#' nrPassThreshold = 1,
#' method = "distclu",
#' maxDist = 20,
#' removeSingletons = TRUE,
#' keepSingletonsAbove = 5,
#' qLow = 0.1,
#' qUp = 0.9,
#' tpmThreshold = 3,
#' plot_lim = c(0, 150),
#' num_core = 4,
#' cager_folder="cager_results/")

cager_clustering <- function(
    ce,
    threshold,
    thresholdIsTpm,
    nrPassThreshold,
    method,
    maxDist,
    removeSingletons,
    keepSingletonsAbove,
    qLow, qUp,
    tpmThreshold,
    plot_lim,
    numCore,
    cagerFolder
){
    multicore <- TRUE
    if(num_core < 2){
        multicore <- FALSE
        num_core <- NULL
    }

    # cluster TSS
    ce <- CAGEr::clusterCTSS(
        ce,
        threshold = 1,
        thresholdIsTpm = TRUE,
        nrPassThreshold = 1,
        method = method,
        maxDist = maxDist,
        removeSingletons = TRUE,
        keepSingletonsAbove = 5,
        useMulticore = multicore,
        nrCores = num_core)

    # calculate IQ range
    ce <- CAGEr::cumulativeCTSSdistribution(
        ce,
        clusters = "tagClusters",
        useMulticore = multicore,
        nrCores = num_core)

    ce <- CAGEr::quantilePositions(
        ce,
        clusters = "tagClusters",
        qLow = qLow,
        qUp = qUp)
    
    # plot IQ width
    pdf(file.path(cagerFolder, "interquartile_width_tagclusters.pdf"))
    my_plot <- plotInterquantileWidth_local(
        ce,
        clusters = "tagClusters",
        tpmThreshold = tpmThreshold,
        qLow = qLow,
        qUp = qUp,
        xlim = plot_lim)
    print(my_plot)
    dev.off()

    # count CTSS
    sample_ctss_count <- list()
    for (sample in CAGEr::sampleLabels(brcage)) {
        sample_ctss_count[[sample]] <- sum(as.vector(CTSStagCountDF(brcage)[[sample]])>0)
    }

    sink(file.path(cagerFolder, "sample_ctss_count.txt"))
    print(data_out)
    sink()

    plot_number_of_ctss(
        sample_ctss_count=sample_ctss_count,
        cagerFolder=cagerFolder,
        yaxistitle="Number of TSS clusters",
        mytitle="CTSS per sample",
        myfilename="ctss_counts_plots.pdf")

    return(ce)
}

