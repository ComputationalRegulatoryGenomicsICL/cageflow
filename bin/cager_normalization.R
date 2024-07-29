#' Normalization
#'
#' @param ce annotated CAGEexp object
#' @param rangeMin range min within which the slope is calculated
#' @param rangeMax range max within which the slope is calculated
#' @param method method of normalizing tag counts
#' @param total_tag_num total number of tags
#' @param cager_folder path to save plots
#' @return ce normalized CAGEexp object
#' @examples
#' cager_normalization(
#' ce,
#' rangeMin=50,
#' rangeMax=100000,
#' method="powerLaw",
#' total_tag_num=1*10^6,
#' cager_folder="cager_results/")

cager_normalization <- function(
    ce,
    rangeMin, rangeMax,
    method,
    total_tag_num,
    cager_folder){

    # merge replicates - automatic extraction of replicate counts
    # only if they follow some pattern
    # ce <- CAGEr::mergeSamples(
    #     ce,
    #     mergeIndex = new_indeces,
    #     mergedSampleLabels = sample_names
    # )

    outlist <- calculateReverseCumulative(
        object = ce,
        values = "raw",
        fitInRange = c(10, 1000))
    tag_count_df <- outlist[[1]]
    slope <- outlist[[2]]
    library_size <- outlist[[3]]
    intercept <- outlist[[4]]

    plots <- plotReverseCumulatives_local(
        tag_count_df=tag_count_df,
        slope=slope,
        intercept=intercept,
        library_size=library_size,
        fitInRange = c(range_min, range_max))
    # TODO: pdf if not too many samples
    png(file.path(cager_folder, "reverse_cumulative.png"))
    print(plots)
    dev.off()

    ce <- normalizeTagCount(
        ce,
        method = method,
        fitInRange = c(range_min, range_max),
        alpha = slope,
        T = total_tag_num)

    return(ce)
}