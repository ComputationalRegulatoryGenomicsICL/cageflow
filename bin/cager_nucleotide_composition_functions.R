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
    bsgenome <- getBSgenome(seq_data)
    ctss_sequences <- lapply(
        ctss_list,
        function(x) {
            BSgenome::getSeq(bsgenome, x)
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
        sample_names) {
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

    return(p)
}

### For dinucleotide composition

expand_ctss_regions <- function(normalized_ctss_list, seq_data) {
    bsgenome <- getBSgenome(seq_data)
    # this might not be necessary
    for (sample in names(normalized_ctss_list)){
        GenomeInfoDb::seqinfo(normalized_ctss_list[[sample]]) <- GenomeInfoDb::seqinfo(bsgenome)[
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
    ctss_dinuc_freq_filt <- lapply(
        ctss_dinuc_freq,
        function(x) {x[names(x) %in% dinuc_values]})
    # convert to dataframe
    cdf_intermediate <- lapply(ctss_dinuc_freq, as.data.frame)
    cdf_almost_good <- dplyr::bind_rows(cdf_intermediate, .id="Name")
    ctss_dinuc_freq_df <- reshape(cdf_almost_good, idvar="Var1", timevar="Name",direction="wide")
    # fill in 0 instead of NA
    ctss_dinuc_freq_df[is.na(ctss_dinuc_freq_df)]<-0
    # order
    ctss_dinuc_freq_df_srt <- arrange(ctss_dinuc_freq_df, desc(across(names(ctss_dinuc_freq_df)[2])))
    # fix names
    names(ctss_dinuc_freq_df_srt) <- gsub("Freq.", "", names(ctss_dinuc_freq_df_srt))
    names(ctss_dinuc_freq_df_srt) <- gsub("Var1", "dinucleotide", names(ctss_dinuc_freq_df_srt))
    rownames(ctss_dinuc_freq_df_srt) <- ctss_dinuc_freq_df_srt$dinucleotide
    return(ctss_dinuc_freq_df_srt)
}

dinucleotide_freq_plot_prep <- function(
    ctss_dinuc_freq_df_tidy){
    # prepare dataframe for ggplot
    ctss_dinuc_freq_df_tidy_gg <- tidyr::pivot_longer(
        ctss_dinuc_freq_df_tidy,
        cols=2:3,
        names_to="samples",
        values_to="percentage"
    )
    # plot initiators as a heatmap - all samples
    # set levels
    ctss_dinuc_freq_df_tidy_gg$dinucleotide <- factor(
        ctss_dinuc_freq_df_tidy_gg$dinucleotide,
        levels=ctss_dinuc_freq_df_tidy$dinucleotide)
    # keep alphabetical order of samples
    column_names <- sort(unique(ctss_dinuc_freq_df_tidy_gg$samples))
    ctss_dinuc_freq_df_tidy_gg$samples <- factor(
        ctss_dinuc_freq_df_tidy_gg$samples,
        levels=column_names)
    
    return(ctss_dinuc_freq_df_tidy_gg)
}

plot_dinucleotide_frequency_heatmap <- function(
        ctss_dinuc_freq_df_tidy_gg) {

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
    return(p)
}

plot_dinucleotide_frequency_histogram <- function(
    ctss_dinuc_freq_df_tidy_gg, col) {

    p <- ggplot(
        data = ctss_dinuc_freq_df_tidy_gg,
        aes(x = dinucleotide, y = percentage, fill = samples)) + 
        scale_fill_manual(values = col) +
        geom_bar(
            stat = "identity",
            position = position_dodge(),
            colour = "black",
            linewidth = 0.25) + 
        coord_flip() +
        xlab("Initiator dinucleotide") + 
        ylab("Percentage") + 
        ggtitle(NULL) + 
        theme_bw() +
        theme(
            text = element_text(size = 14, colour = "black"), 
            axis.text.x = element_text(size = 14, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        scale_y_continuous(limits = c(0, 30)) +
        labs(fill = NULL)

    return(p)
}

plot_dinucleotide_frequency <- function(
        ctss_dinuc_freq_df_tidy) {
    
    ctss_dinuc_freq_df_tidy_gg <- dinucleotide_freq_plot_prep(ctss_dinuc_freq_df_tidy)

    column_names <- sort(unique(ctss_dinuc_freq_df_tidy_gg$samples))
    col = viridis::magma(
        length(column_names),
        alpha = 0.8)[length(column_names):1]

    # heatmap if more than 10 samples, otherwise barplot
    if (length(column_names) > 10){
        p <- plot_dinucleotide_frequency_heatmap(
            ctss_dinuc_freq_df_tidy_gg=ctss_dinuc_freq_df_tidy_gg
        )
    } else {
        p <- plot_dinucleotide_frequency_histogram(
            ctss_dinuc_freq_df_tidy_gg=ctss_dinuc_freq_df_tidy_gg,
            col=col)
    }
    return(p)
}
