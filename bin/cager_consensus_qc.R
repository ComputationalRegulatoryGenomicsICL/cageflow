#' QC plots on consensus CTSS across samples
#'
#' @param ce normalized CAGEexp object
#' @param pcarank Rank of PCA
#' @examples
#' consensus_qc(
#' ce,
#' pcarank = 50)

consensus_qc <- function(
        ce,
        pcarank){

    consclustTmp <- CAGEr::consensusClustersTpm(ce)
    # save matrix of sample per consensus cluster TPM
    write.table(
        consclustTmp,
        file="consensus_clusters_tpm.csv")
    print("Consensus cluster tpms saved")

    # count and plot the number of consensus clusters with signal
    sample_cons_ctss_count <- list()
    for (sample in CAGEr::sampleLabels(ce)) {
        sample_cons_ctss_count[[sample]] <- sum(
            as.vector(consclustTmp[,sample]) > 0)
    }
    sample_cons_ctss_count[["Union"]] <- dim(consclustTmp)[1]
    consensus_ctss_plot <- plot_number_of_ctss(
        sample_ctss_count=sample_cons_ctss_count,
        yaxistitle="Number of non-zero consensus clusters",
        mytitle="Non-zero consensus clusters per sample and union")
    save_plot(
        "consensus_counts_plot.pdf",
        consensus_ctss_plot)
    print("Number of non-zero consensus clusters plotted")

    # plot PCAs
    pca_plot <- plot_pcs(ce_tmp=consclustTmp, pcarank=pcarank)
    save_plot(
        "consensus_clusters_pca_plot.pdf",
        pca_plot)
    print("Consensus cluster PCA plotted")
}