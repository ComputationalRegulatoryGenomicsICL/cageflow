# Functions to calulcate and plot dinucleotide composition

extract_dinucleotide_information <- function(ce, reference_name) {
    bsgenome <- BSgenome::getBSgenome(reference_name)
    sample_names <- CAGEr::sampleLabels(ce)
    weigthed_dinuc_vals <- list()
    for (sample in sample_names) {
        tmp <- as.data.frame(tagClustersGR(
            ce,
            sample = sample,
            qLow = 0.1,qUp = 0.9))
        tmp <- GRanges(
            seqnames = tmp$seqnames,
            ranges =  IRanges(
                start = tmp$dominant_ctss.pos,
                end = tmp$dominant_ctss.pos),
            strand = tmp$strand,
            score=tmp$dominant_ctss.score,
            seqlengths = seqlengths(bsgenome))
        tmp <- IRanges::promoters(
            tmp,
            upstream = 1,
            downstream = 1)
        tmp <- tmp[width(trim(tmp)) == 2]
        tmp$dinucleotide <- as.data.frame(getSeq(bsg,tmp))$x

        dinuc_n_score<-as.data.frame(tmp)[,c("score","dinucleotide")]
        dinucleotide<-unique(dinuc_n_score$dinucleotide)
        dinuc_vals<-data.frame(
            dinucleotide=dinucleotide,
            sum_score=NA)
        for (l in dinucleotide) {
            
            dinuc_vals[which(dinuc_vals$dinucleotide==l),]$sum_score <- sum(
                dinuc_n_score[which(dinuc_n_score$dinucleotide==l),]$score)
        
        }
        
        dinuc_vals$proportion <- dinuc_vals$sum_score/sum(dinuc_vals$sum_score)
        
        dinuc_vals<-dinuc_vals[
            order(dinuc_vals$proportion, decreasing=F),]
        dinuc_vals$dinucleotide <- factor(
            dinuc_vals$dinucleotide,
            levels = dinuc_vals$dinucleotide)

        weigthed_dinuc_vals[[sample]] <- dinuc_vals
    }
    return(weigthed_dinuc_vals)
}

# count_dinucleotide_frequency <- function(ctss_sequences) {
#     # count dinucleotides in initiator CTSS
#     ctss_dinuc_count <- lapply(ctss_sequences, table)
#     # convert dinucleotide count to percentages
#     ctss_dinuc_freq <- lapply(
#         ctss_dinuc_count,
#         function(x) {x/sum(x) * 100}
#     )
#     # set levels of dinucleotide usage for plotting
#     all_dinucleotides <- unique(unlist(lapply(ctss_dinuc_freq, names)))
#     # remove dinuc with N base
#     all_nonn_dinuc <- all_dinucleotides[!grepl("N", all_dinucleotides)]
#     # remove single letter long "dinucleotides"
#     true_dinuc <- unlist(lapply(all_nonn_dinuc, function(x) {x[nchar(x) > 1]}))
#     dinuc_values <- true_dinuc[order(true_dinuc)]
#     ctss_dinuc_freq_filt <- lapply(
#         ctss_dinuc_freq,
#         function(x) {x[names(x) %in% dinuc_values]})
#     # convert to dataframe
#     cdf_intermediate <- lapply(ctss_dinuc_freq, as.data.frame)
#     cdf_almost_good <- dplyr::bind_rows(cdf_intermediate, .id="Name")
#     ctss_dinuc_freq_df <- reshape(cdf_almost_good, idvar="Var1", timevar="Name",direction="wide")
#     # fill in 0 instead of NA
#     ctss_dinuc_freq_df[is.na(ctss_dinuc_freq_df)]<-0
#     # order
#     ctss_dinuc_freq_df_srt <- arrange(ctss_dinuc_freq_df, desc(across(names(ctss_dinuc_freq_df)[2])))
#     # fix names
#     names(ctss_dinuc_freq_df_srt) <- gsub("Freq.", "", names(ctss_dinuc_freq_df_srt))
#     names(ctss_dinuc_freq_df_srt) <- gsub("Var1", "dinucleotide", names(ctss_dinuc_freq_df_srt))
#     rownames(ctss_dinuc_freq_df_srt) <- ctss_dinuc_freq_df_srt$dinucleotide
#     return(ctss_dinuc_freq_df_srt)
# }

# dinucleotide_freq_plot_prep <- function(
#     ctss_dinuc_freq_df_tidy){
#     # prepare dataframe for ggplot
#     ctss_dinuc_freq_df_tidy_gg <- tidyr::pivot_longer(
#         ctss_dinuc_freq_df_tidy,
#         cols=2:length(ctss_dinuc_freq_df_tidy),
#         names_to="samples",
#         values_to="percentage"
#     )
#     # plot initiators as a heatmap - all samples
#     # set levels
#     ctss_dinuc_freq_df_tidy_gg$dinucleotide <- factor(
#         ctss_dinuc_freq_df_tidy_gg$dinucleotide,
#         levels=ctss_dinuc_freq_df_tidy$dinucleotide)
#     # keep alphabetical order of samples
#     column_names <- sort(unique(ctss_dinuc_freq_df_tidy_gg$samples))
#     ctss_dinuc_freq_df_tidy_gg$samples <- factor(
#         ctss_dinuc_freq_df_tidy_gg$samples,
#         levels=column_names)
    
#     return(ctss_dinuc_freq_df_tidy_gg)
# }

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
