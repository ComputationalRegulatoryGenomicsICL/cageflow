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

# Helper function to Pre-calculate a vector of correlation coefficients
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

# Re-implement the pairs function to prevent coercion to data.frame
pairs.DataFrame <- function (x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc, 
    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3, 
    cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1, 
    log = "") 
{
    if (doText <- missing(text.panel) || is.function(text.panel)) 
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
            y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
        oma, ...) {
        xpd <- NA
        if (side%%2L == 1L && xl[j]) 
            xpd <- FALSE
        if (side%%2L == 0L && yl[i]) 
            xpd <- FALSE
        if (side%%2L == 1L) 
            Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    # if (!is.matrix(x)) {
    #     x <- as.data.frame(x)
    #     for (i in seq_along(names(x))) {
    #         if (is.factor(x[[i]]) || is.logical(x[[i]])) 
    #             x[[i]] <- as.numeric(x[[i]])
    #         if (!is.numeric(unclass(x[[i]]))) 
    #             stop("non-numeric argument to 'pairs'")
    #     }
    # }
    # else if (!is.numeric(x)) 
    #     stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2L) 
        stop("only one column in the argument to 'pairs'")
    if (!all(horInd >= 1L & horInd <= nc)) 
        stop("invalid argument 'horInd'")
    if (!all(verInd >= 1L & verInd <= nc)) 
        stop("invalid argument 'verInd'")
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) 
                labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels)) 
            doText <- FALSE
    }
    oma <- if ("oma" %in% nmdots) 
        dots$oma
    main <- if ("main" %in% nmdots) 
        dots$main
    if (is.null(oma)) 
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
        4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log)) 
        xl[log] <- yl[log] <- TRUE
    else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    for (i in if (row1attop) 
        verInd
    else rev(verInd)) for (j in horInd) {
        l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
            ""))
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ..., log = l)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
                localAxis(1L + 2L * row1attop, x[, j], x[, i], 
                  ...)
            if (i == nc && (j%%2L || !has.upper || !has.lower)) 
                localAxis(3L - 2L * row1attop, x[, j], x[, i], 
                  ...)
            if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
                localAxis(2L, x[, j], x[, i], ...)
            if (j == nc && (i%%2L || !has.upper || !has.lower)) 
                localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  localDiagPanel(x[, i], ...)
                if (doText) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  xlp <- if (xl[i]) 
                    10^0.5
                  else 0.5
                  ylp <- if (yl[j]) 
                    10^label.pos
                  else label.pos
                  text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) 
                localLowerPanel(x[, j], x[, i], ...)
            else localUpperPanel(x[, j], x[, i], ...)
            if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
            font = font.main)
    }
    invisible(NULL)
}


plotCorrelation2_local <- function( expr.table, samples, method
                               , tagCountThreshold, applyThresholdBoth
                               , digits, plot_pairs=FALSE) {
  # Select samples
  if(samples == "all"){
    samples <- colnames(expr.table)
  } else if (all(samples %in% colnames(expr.table))) {
    expr.table <- expr.table[,samples]
  } else stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
  nr.samples <- length(samples)
  
  # Pre-calculate a vector of correlation coefficients
  corr.v <- corVector(expr.table, method, tagCountThreshold, applyThresholdBoth)
  
  # Add pseudocount to null values so that the plot axes are correctly set.
  pseudocount <- min(sapply(expr.table, function(x) min(x[x>0]))) / 2
  expr.table  <- DataFrame(lapply( expr.table
                                   , function(x) {x[x==0] <- pseudocount ; x}))
  
  # This closure retreives correlation coefficients one after the other.
  mkPanelCor <- function() {
    i <- 1
    function(x, y, digits=2, prefix="", cex.cor, ...) {
      r <- corr.v[i]
      i <<- i + 1
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      gmean <- memoise(function(x) {
        return(exp(mean(log(c(pseudocount, max(x))))))
      })
      text(gmean(x), gmean(y), txt, cex = cex.cor * sqrt(r))
    }
  }
  panel.cor <- mkPanelCor()
  
  # uniqueSignif returns a plain data.frame, because compression is already
  # maximal in this context.  See benchmarking of alternatives algorithms
  # in the file "benchmarks/unique-signif.md" in the CAGEr's Git repository.
  uniqueSignif <- function(x, y, digits = 0, log = c("", "xy")) {
    log <- match.arg(log)
    if (log == "xy") {x <- log1p(x) ; y <- log1p(y)}
    u  <- unique(Rle(complex( real      = decode(signif(x, digits = digits))
                              , imaginary = decode(signif(y, digits = digits)))))
    df <- data.frame(x = Re(u), y = Im(u))
    if (log == "xy") df <- expm1(df)
    df
  }
  
  # Thresholds are lowered of a minute amont because the rounding in the log
  # scale in uniqueSignif adds minute errors to values that were already round.
  pointsUnique <- function(x,y,...) {
    df <- uniqueSignif(x, y, digits = digits, log = "xy")
    df <- .applyThreshold(df, tagCountThreshold * 0.999, applyThresholdBoth)
    df <- df[rowSums(df) > pseudocount * 1.999,] # Remove the (0,0) point.
    points(df, ...)
    # p <- ggplot(df, aes(x=x, y=y)) + 
    #   geom_density_2d(aes(fill = ..level..), geom = "polygon") + 
    #   theme_classic()
    # print(p)
  }
  
  if (plot_pairs){
    pdf("plots/correlations_plot.pdf")
    pairs( expr.table
           , lower.panel = pointsUnique
           , upper.panel = panel.cor
           , pch = "."
           , cex = 4
           , log = "xy"
           , las = 1
           , xaxp = c(1,10,1)
           , yaxp = c(1,10,1)
           , labels = samples)
    dev.off()
    pairs_plot <- list(expr.table, pointsUnique, panel.cor, samples, corr.v, pseudocount)
    saveRDS(pairs_plot, "plots/correlations_plot.rds")
  }

  p <- ggpairs(
    expr.table,
    upper = list(continuous = "cor"),
    lower = list(continuous = "density"),
    diag = list(continuous = "densityDiag")
  )
  ggsave("ggpairs_test.pdf", p)
  
  # Return a correlation matrix
  corr.m <- matrix(1, nr.samples, nr.samples)
  colnames(corr.m) <- samples
  rownames(corr.m) <- samples
  corr.m[lower.tri(corr.m)] <- corr.v
  corr.m[upper.tri(corr.m)] <- t(corr.m)[upper.tri(corr.m)]
  return(corr.m)
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
