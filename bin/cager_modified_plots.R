plotInterquantileWidth_local <- function(object, clusters, tpmThreshold, qLow, qUp, xlim){
#  clusters <- match.arg(clusters)
#	getClustFun <- switch( clusters
#	                     , tagClusters       = tagClustersGR
#	                     , consensusClusters = consensusClustersGR)
  getClustFun <- tagClustersGR 
  # Extract a list of data frames in "long" format for ggplot
	iqwidths <- lapply(seq_along(sampleLabels(object)), function(x) {
    gr <- getClustFun(object, x, qLow = qLow, qUp = qUp)
    gr <- gr[score(gr) >= tpmThreshold]
    data.frame(
      sampleName  = sampleLabels(object)[[x]],
      iq_width    = decode(gr$interquantile_width))
	})
	
	# Bind them together and set factor labels in proper order
	iqwidths <- do.call(rbind, iqwidths)
  iqwidths$sampleName <- factor(iqwidths$sampleName, levels = sampleLabels(object))
  iqwidths <- iqwidths[iqwidths$iq_width >= xlim[1] & iqwidths$iq_width <= xlim[2],]
  
	# binsize <- round(max(iqwidths$iq_width)/2)
  binsize <- ceiling(max(iqwidths$iq_width)/2)

  # debug
  # cat("nrow(iqwidths) = ", nrow(iqwidths), "\n")
  # cat("max(iqwidths$iq_width) = ", max(iqwidths$iq_width), "\n")
  # cat("binsize = ", binsize, "\n")
  # end_debug
	
	iqwidth_plot <- ggplot2::ggplot(iqwidths) +
	  ggplot2::aes_string(x = "iq_width") +
	  ggplot2::scale_fill_manual(values = names(sampleLabels(object))) +
	  ggplot2::geom_histogram(bins = binsize) +
	  ggplot2::facet_wrap("~sampleName") +
	  ggplot2::ggtitle(paste0(
	    switch(clusters, tagClusters = "Tag Clusters", consensusClusters = "Consensus Clusters"),
	    " interquantile width (quantile ", qLow, " to ", qUp, ")")) +
	  ggplot2::xlab("Interquantile width (bp)") +
	  ggplot2::ylab("Frequency")

	return(iqwidth_plot)
}

# Helper function to apply threshold pairwise
.applyThreshold <- function(df, tagCountThreshold, applyThresholdBoth) {
  if (applyThresholdBoth) {
    idx <- (df[[1]] >= tagCountThreshold) & (df[[2]] >= tagCountThreshold)
  } else {
    idx <- (df[[1]] >= tagCountThreshold) | (df[[2]] >= tagCountThreshold)
  }
  df[idx,]
}

# # Helper function to Pre-calculate a vector of correlation coefficients
corVector <- function(expr.table, method, tagCountThreshold, applyThresholdBoth) {
  corTreshold <- function(x, y, method) {
    df <- data.frame(x, y)
    df <- .applyThreshold(df, tagCountThreshold, applyThresholdBoth)
    cor(x = df$x, y = df$y, method = method)
  }
  nr.samples <- ncol(expr.table)
  corr.v <- numeric()
  for (i in 1:(nr.samples-1)) {
    for (j in (min(i+1, nr.samples)):nr.samples) {
      corr.v <- append(corr.v, corTreshold(expr.table[[i]], expr.table[[j]], method))
    }
  }
  corr.v
}


calculate_correlation_matrix <- function(
    expr.table, samples, method, tagCountThreshold,
    applyThresholdBoth){

  # Select samples
  if(samples == "all"){
    samples <- colnames(expr.table)
  } else if (all(samples %in% colnames(expr.table))) {
    expr.table <- expr.table[,samples]
  } else stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
  nr.samples <- length(samples)

  # Pre-calculate a vector of correlation coefficients
  corr.v <- corVector(expr.table, method, tagCountThreshold, applyThresholdBoth)

  corr.m <- matrix(1, nr.samples, nr.samples)
  colnames(corr.m) <- samples
  rownames(corr.m) <- samples
  corr.m[lower.tri(corr.m)] <- corr.v
  corr.m[upper.tri(corr.m)] <- t(corr.m)[upper.tri(corr.m)]

}


.fit.power.law.to.reverse.cumulative <- function(values, val.range = c(10, 1000)) {
  
  # using data.table package	
  v <- data.table(num = 1, nr_tags = values)
  num <- nr_tags <- NULL # Keep R CMD check happy.
  v <- v[, sum(num), by = nr_tags]
  setkeyv(v, "nr_tags")
  v$V1 <- rev(cumsum(rev(v$V1)))
  setnames(v, c('nr_tags', 'reverse_cumulative'))
  
  # check if range values are > 0
  if(any(val.range < 1)){
    stop(paste("The range of values for fitting the power law",
               "arg 'fitInRange', expects integers > 0."))
  }
  
  v <- v[nr_tags >= min(val.range) & nr_tags <= max(val.range)]
  
  # check if specified range values have at least 1 entry in v 
  if(nrow(v) < 1){
    stop(paste("Selected range for fitting the power law does not contain any",
               "tag count values. Consider changing/increasing 'fitInRange'."))
  }
  #	v <- aggregate(values, by = list(values), FUN = length)
  #	v$x <- rev(cumsum(rev(v$x)))
  #	colnames(v) <- c('nr_tags', 'reverse_cumulative')
  #	v <- subset(v, nr_tags >= min(val.range) & nr_tags <= max(val.range))
  
  lin.m <- lm(log(reverse_cumulative) ~ log(nr_tags), data = v)
  
  a <- coefficients(lin.m)[2]
  b <- coefficients(lin.m)[1]
  
  # check if specified range values have >1 entries in v 
  if(is.na(a) && b == 0){
    stop(paste("Selected range for fitting the power law does not", 
               "contain enough values. Consider changing/increasing 'fitInRange'."))
  }
  return(c(a, b))
}


prepare_counts <- function(tag_counts){
  tag_xy <- tibble(name=character(), x=numeric(), y=numeric())
  # for each column, add x and y values
  for (column_idx in 1:dim(tag_counts)[2]){
    tag_column <- tag_counts[,column_idx]
    values <- sort(Rle(tag_column), decreasing = TRUE)
    values <- values[values != 0]
    
    x_values <- runValue(values)
    y_values <- cumsum(runLength(values))
    tag_xy <- tag_xy %>% add_row(
      name= colnames(tag_counts)[column_idx], x=x_values, y=y_values
    )
  }
  return(tag_xy)
}

calculateReverseCumulative <- function(
    object,
    values = c("raw", "normalized"),
    fitInRange = c(10, 1000)
    ){
    sample.labels <- sampleLabels(object)
    values <- match.arg(values)
    old.par <- par(mar = c(5,5,5,2))
    on.exit(par(old.par))
    cols <- names(sample.labels)
    
    tag.count <- switch( values
                        , raw        = CTSStagCountDF(object)
                        , normalized = CTSSnormalizedTpmDF(object))
    
    if(! is.null(fitInRange)) {
        fit.coefs.m <- as.matrix(data.frame(lapply(tag.count, function(x) {
        .fit.power.law.to.reverse.cumulative(values = decode(x), val.range = fitInRange)})))
        fit.slopes <- fit.coefs.m[1,]
        names(fit.slopes) <- sample.labels
        reference.slope <- min(median(fit.slopes), -1.05)
        reference.library.size <- 10^floor(log10(median(sapply(tag.count, sum))))
        #reference.intercept <- log(reference.library.size/zeta(-1*reference.slope))  # intercept on natural logarithm scale
        reference.intercept <- log10(reference.library.size/VGAM::zeta(-1*reference.slope))  # intercept on log10 scale used for plotting with abline
    }
    tag_count_df <- prepare_counts(tag.count)

    # return(list(tag_count_df, reference.slope, reference.library.size, reference.intercept))
    return(list(tag_count_df, reference.slope, reference.library.size, reference.intercept, fit.slopes))

}

plotReverseCumulatives_local <- function(
    tag_count_df,
    slope, 
    intercept,
    library_size,
    fit.slopes,
    fitInRange = c(10, 1000),
    main = NULL, legend = TRUE,
    xlab = "number of CAGE tags", ylab = "number of CTSSs (>= nr tags)",
    xlim = c(1, 1e5), ylim = c(1, 1e6)){
  
  # plot nicely
  # facets grid
  # remove alpha, but add as title for alpha and subtitle for T
  # geom_vline for vertical as two numbers
  # diagonal as separate geom_abline slope and intercept
  
  plot_out <- ggplot2::ggplot(tag_count_df) +
    ggplot2::aes(x=x, y=y) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(. ~name) +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(title="Reference distribution:",
         subtitle = paste0("alpha= ", sprintf("%.2f", -1*slope), " T= ", library_size), # library.size
         x =xlab, y = ylab) +
    ggplot2::geom_text(data = tag_count_df,
       mapping= aes(
         x=10, y= 10,
         label = paste0("alpha= ", formatC(-1*fit.slopes[name], format = "f", digits = 2)))) +
    ggplot2::geom_vline(xintercept=fitInRange, linetype="dotted") +
    ggplot2::geom_abline(slope = slope, intercept = intercept,
                         linetype="longdash", colour="#7F7F7F7F")
    
  return(plot_out)
}
