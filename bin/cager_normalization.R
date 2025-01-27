#' Normalization
#'
#' @param ce annotated CAGEexp object
#' @param rangeMin range min within which the slope is calculated
#' @param rangeMax range max within which the slope is calculated
#' @param method method of normalizing tag counts
#' @param T_norm total number of tags
#' @return ce normalized CAGEexp object
#' @examples
#' cager_normalization(
#' ce,
#' rangeMin=10,
#' rangeMax=10000,
#' method="powerLaw",
#' T_norm=1*10^6)

cager_normalization <- function(
    ce,
    rangeMin, rangeMax,
    method,
    T_norm){

    outlist <- calculateReverseCumulative(
        object = ce,
        values = "raw",
        fitInRange = c(rangeMin, rangeMax))
    tag_count_df <- outlist[[1]]
    slope <- outlist[[2]]
    library_size <- outlist[[3]]
    intercept <- outlist[[4]]
    fit.slopes <- outlist[[5]]

    revcum_plots <- plotReverseCumulatives_local(
        tag_count_df=tag_count_df,
        slope=slope,
        intercept=intercept,
        library_size=library_size,
        fit.slopes = fit.slopes,
        fitInRange = c(rangeMin, rangeMax))
    
    save_plot(
        "reverse_cumulative_plot.pdf",
        revcum_plots)

    ce <- CAGEr::normalizeTagCount(
        ce,
        method = method,
        fitInRange = c(rangeMin, rangeMax),
        alpha = -slope,
        T = T_norm)

    return(ce)
}