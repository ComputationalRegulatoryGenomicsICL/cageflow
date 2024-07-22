#' Initial quality control of CAGE reads
#'
#' @param ce CAGEexp object with tags from BAM or BigWig
#' @param tx_annotation Annotation for TxDb object, eg TxDb.Hsapiens.UCSC.hg38.knownGene
#' @param cager_folder Path to save the plots 
#' @examples
#' cager_qc(
#' ce,
#' tx_annotation=TxDb.Hsapiens.UCSC.hg38.knownGene,
#' cager_folder=cager_results/ )

cager_qc <- function(
    ce,
    tx_annotation,
    cager_folder){

    pdf(file.path(cager_folder, "tag_region_annotation.pdf"))
    annotations <- CAGEr::plotAnnot(ce, "counts")
    print(annotations)
    dev.off()

    corr_m <- plotCorrelation2_local(
        CTSStagCountDF(ce),
        samples = "all",
        tagCountThreshold = 1,
        applyThresholdBoth = FALSE,
        method = "pearson",
        digits = 3,
        toPlot = FALSE)
    
    # save intermediate file
    saveRDS(corr_m, file.path(cager_folder, "corr_m.rds"))

    # plot correlations in heatmap format
    pdf(file.path(cager_folder, "correlations_heatmap.pdf"))
    heatmap.2(corr_m, trace="none", margins=c(12, 12),cexRow=0.2)
    dev.off()

}