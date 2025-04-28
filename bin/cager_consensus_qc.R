#' QC plots on consensus CTSS across samples
#'
#' @param ce normalized CAGEexp object
#' @param heatmap_cex Text size for heatmap
#' @examples
#' consensus_qc(
#' ce,
#' heatmap_cex=12)

consensus_qc <- function(
        ce,
        heatmap_cex){

    consclustTpm <- CAGEr::consensusClustersTpm(ce)
    # save matrix of sample per consensus cluster TPM
    write.table(
        consclustTpm,
        file="tables/consensus_clusters_tpm.csv")
    print("Consensus cluster tpms saved")

    # count and plot the number of consensus clusters with signal
    sample_cons_ctss_count <- list()
    for (sample in CAGEr::sampleLabels(ce)) {
        sample_cons_ctss_count[[sample]] <- sum(
            as.vector(consclustTpm[,sample]) > 0)
    }
    sample_cons_ctss_count[["Union"]] <- dim(consclustTpm)[1]
    consensus_ctss_plot <- plot_number_of_tag_clusters(
        sample_tag_count=sample_cons_ctss_count,
        yaxistitle="Number of non-zero consensus clusters",
        mytitle="Non-zero consensus clusters per sample and union")
    save_plot(
        "consensus_counts_plot.pdf",
        consensus_ctss_plot)
    print("Number of non-zero consensus clusters plotted")

    # plot correlation
    plot_correlation(
        datatype="consensus_clusters",
        dataframe=consclustTpm,
        corrplot_tagCountThreshold=0,
        heatmap_cex=heatmap_cex)
    print("CTSS correlation plotted")

    # plot PCAs
    pca_plot <- plot_pcs(consclustTpm)
    save_plot(
        "consensus_clusters_pca_plot.pdf",
        pca_plot)
    print("Consensus cluster PCA plotted")
}