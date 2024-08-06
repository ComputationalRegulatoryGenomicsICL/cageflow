# Functions to calulcate and plot nucleotide and dinucleotide compositions

extract_ctss_normalized_tmp_per_sample <- function(ce, tpm_threshold = 1) {
    sample_names <- sampleLabels(ce)
    normalized_ctss_list <- list()
    for (sample in sample_names) {
        unfiltered <- CAGEr::CTSSnormalizedTpmGR(ce, sample)
        normalized_ctss_list[[sample]] <- unfiltered[score(unfiltered) >= tpm_threshold]
    }
    return(normalized_ctss_list)
}

extract_ctss_sequences <- function(ctss_list, seq_data) {
    ctss_sequences <- lapply(
        ctss_list,
        function(x) {
            BSgenome::getSeq(seq_data, x)
        }
    )
    return(ctss_sequences)
}

### For nucleotide composition

calculate_nucleotide_frequency <- function(ctss_sequences) {
    sample_names <- names(ctss_sequences)
    ctss_nucl_freq <- lapply(
        ctss_sequences,
        function(x) {
            Biostrings::nucleotideFrequencyAt(x, at = 1) / length(x) * 100
        }
    )
    ctss_nucl_freq_df <- data.frame(
        matrix(
            unlist(ctss_nucl_freq),
            nrow = length(ctss_nucl_freq),
            byrow = TRUE
        )
    )
    colnames(ctss_nucl_freq_df) <- names(ctss_nucl_freq[[1]])
    rownames(ctss_nucl_freq_df) <- sample_names
    ctss_nucl_freq_df$samples <- rownames(ctss_nucl_freq_df)
    ctss_nucl_freq_df_tidy <- tidyr::gather(
        ctss_nucl_freq_df,
        Feature,
        Percentage,
        (c("A", "C", "G", "T"))
    )
    return(list(ctss_nucl_freq_df_tidy, sample_names))
}


plot_nucleotide_frequency <- function(
        ctss_nucl_freq_df_tidy,
        sample_names,
        outfilepath,
        pdfheight=3,
        pdfwidth=5) {
    # plot start nucleotide frequencies as a histogram
    # set levels, alphabetical order of samples
    ctss_nucl_freq_df_tidy$samples <- factor(
        ctss_nucl_freq_df_tidy$samples,
        levels = sample_names)

    # set levels for nucleotides
    # reverse order to plot nicer
    ctss_nucl_freq_df_tidy$Feature <- factor(
        ctss_nucl_freq_df_tidy$Feature,
        levels = c("T", "G", "C", "A"))

    # plot start nucleotide frequencies as a histogram for all samples
    # select colours
    col <- viridis::magma(5, alpha = 0.7)[c(1,4,3,2)]

    # select colours - standard nucleotide color scheme 
    # A, C, G, T (green, blue, yellow, red)
    col <- c("darkolivegreen3", "cornflowerblue", "gold", "firebrick3")[4:1]

    p <- ggplot(
            data = ctss_nucl_freq_df_tidy,
            aes(x = samples, y = Percentage, fill = Feature)) +
        geom_bar(stat = "identity", width = 0.75, colour = "black", linewidth = 0.125) +
        coord_flip() +
        scale_fill_manual("Nucleotide", values = col) +
        theme_bw() +
        theme(text = element_text(size = 14, colour = "black"),
            legend.title = element_blank(),
            axis.title.x = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        labs(y = "Percentage", x = NULL)

    pdf(file=outfilepath, height = pdfheight, width = pdfwidth)
    print(p)
    dev.off()
    return("Nucleotide frequencies plotted")
}

### For dinucleotide composition

expand_ctss_regions <- function(normalized_ctss_list) {
    # this might not be necessary
    for (sample in names(normalized_ctss_list)){
        GenomeInfoDb::seqinfo(normalized_ctss_list[[sample]]) <- GenomeInfoDb::seqinfo(seq_data)[
            GenomeInfoDb::seqlevels(normalized_ctss_list[[sample]])
        ]
    }
    # expand to include upstream nucleotide (to get (-1, 0))
    # here filtering for promoters
    # From IRanges documentation:
    # The full range is defined as, (start(x) -upstream) to (start(x) + downstream - 1)
    expanded_ctss_list <- lapply(
        normalized_ctss_list,
        function(x) {
            IRanges::promoters(x, upstream = 1, downstream = 1)
        }
    )
    for (sample in names(expanded_ctss_list)) {
        expanded_ctss_list[[sample]] <- GenomicRanges::trim(
            expanded_ctss_list[[sample]])
    }
    
    return(expanded_ctss_list)
}

count_dinucleotide_frequency <- function(ctss_sequences) {
    # count dinucleotides in initiator CTSS
    ctss_dinuc_count <- lapply(ctss_sequences, table)
    # convert dinucleotide count to percentages
    ctss_dinuc_freq <- lapply(
        ctss_dinuc_count,
        function(x) {x/sum(x) * 100}
    )
    # set levels of dinucleotide usage for plotting
    all_dinucleotides <- unique(unlist(lapply(ctss_dinuc_freq, names)))
    # remove dinuc with N base
    all_nonn_dinuc <- all_dinucleotides[!grepl("N", all_dinucleotides)]
    # remove single letter long "dinucleotides"
    true_dinuc <- unlist(lapply(all_nonn_dinuc, function(x) {x[nchar(x) > 1]}))
    dinuc_values <- true_dinuc[order(true_dinuc)]
    # ensure the lists are of the same length - remove entries with N base or single base
    ctss_dinuc_freq_filt <- lapply(
        ctss_dinuc_freq,
        function(x) {x[names(x) %in% dinuc_values]})
    # set levels in value order based on first entry order
    dinuc_levels <- names(sort(ctss_dinuc_freq_filt[[1]]))
    # convert to dataframe
    ctss_dinuc_freq_df <- as.data.frame(ctss_dinuc_freq_filt)
    # tidy the dataframe
    ctss_dinuc_freq_df_tidy <- ctss_dinuc_freq_df[, grep(
        x = colnames(ctss_dinuc_freq_df),
        pattern = "Freq")]
    # attach dinucleotide information
    ctss_dinuc_freq_df_tidy <- cbind(
        ctss_dinuc_freq_df_tidy,
        ctss_dinuc_freq_df[,1])
    # rename columns
    column_names <- names(ctss_sequences)
    colnames(ctss_dinuc_freq_df_tidy) <- c(column_names, "dinucleotide")
    rownames(ctss_dinuc_freq_df_tidy) <- ctss_dinuc_freq_df_tidy$dinucleotide
    return(list(ctss_dinuc_freq_df_tidy, dinuc_levels, column_names))
}

plot_dinucleotide_frequency <- function(
        ctss_dinuc_freq_df_tidy,
        dinuc_levels,
        column_names,
        outfilepath,
        pdfheight = 12,
        pdfwidth = 10) {
    # prepare dataframe for ggplot
    ctss_dinuc_freq_df_tidy_gg <- tidyr::gather(
        ctss_dinuc_freq_df_tidy,
        samples,
        percentage,
        colnames(ctss_dinuc_freq_df_tidy[,1:(length(ctss_dinuc_freq_df_tidy)-1)]))

    # plot initiators as a histogram - all samples
    # set levels
    ctss_dinuc_freq_df_tidy_gg$dinucleotide <- factor(
        ctss_dinuc_freq_df_tidy_gg$dinucleotide,
        levels = dinuc_levels)
    # keep alphabetical order of samples
    ctss_dinuc_freq_df_tidy_gg$samples <- factor(
        ctss_dinuc_freq_df_tidy_gg$samples,
        levels = column_names)

    col = viridis::magma(
        length(column_names),
        alpha = 0.8)[length(column_names):1]

    p <- ggplot(
        data = ctss_dinuc_freq_df_tidy_gg,
        aes(x = dinucleotide, y = samples, fill = percentage)) + 
        geom_tile() + 
        xlab("Initiator dinucleotide") + 
        ylab("Samples") + 
        ggtitle(NULL) + 
        theme_bw() +
        theme(
            text = element_text(size = 4, colour = "black"), 
            axis.text.x = element_text(size = 4, colour = "black",angle = 90),
            axis.text.y = element_text(size = 3, colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        labs(fill = "%")
    # todo: fix this with inputs
    # pdfheight <- (length(colnames(ctss_dinuc_freq_df_tidy))-1)*0.4+2
    pdf(file = outfilepath, height = pdfheight, width = pdfwidth)
    print(p)
    invisible(dev.off())
    print("Dinucleotide frequency plotted")
}


