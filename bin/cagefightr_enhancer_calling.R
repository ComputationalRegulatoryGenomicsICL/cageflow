#' Code of exporting CAGEexp to dgCMatrix is taken from CAGEr 2.12.0
#' Enhancer calling
#'
#' @param ce initial CAGEexp object with CTSS values
#' @param cfBalanceThreshold threshold for the cagefightr balance score
#' @param txdb TxDb object for annotation
#' @return enhancers
#' @examples
#' cagefightr_enhancers(
#' ce,
#' cfBalanceThreshold = 0.95,
#' txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
#' )

cagefightr_enhancers <- function(
    ce,
    cfBalanceThreshold,
    txdb){

    se <- CTSStagCountSE(ce) # TODO: take normalizedSE
    colData(se) <- colData(ce)
    rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
    colData(se)$Name <- colData(se)$sampleLabels
    assays(se, withDimnames=FALSE) <- List(counts=as(as.matrix(as.data.frame(assay(se))), "dgCMatrix"))
    # . 
    # .
    # . from Damir
    BCs <- clusterBidirectionally(se, balanceThreshold=cfBalanceThreshold)
    # Calculate number of bidirectional samples
    BCs <- calcBidirectionality(BCs, samples=se)
    # remove BCs not observed to be bidirectional in one or more samples
    enhancers <- subset(BCs, bidirectionality > 0)
    # annotate enhancers and remove those overlapping known transcripts
    enhancers <- assignTxType(enhancers, txModels=txdb)
    enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))
    return(enhancers)
}


identify_sample_specific_enhancers <- function(BCs){
    
    supported_enhancers <- CAGEfightR::subsetBySupport(
        BCs,
        inputAssay="counts",
        unexpressed=0,
        minSamples=1)
    sample_per_enhancer <- supported_enhancers@assays@data$counts
    return(sample_per_enhancer)
}

save_each_sample_to_bed <- function(sample_per_enhancer, cagerFolder){
    # this function also calculates enhancer widths
    enhancer_widths <- list()
    for (sample_name in colnames(sample_per_enhancer)){
        enhancer_widths[[sample_name]] <- c()
        sample_bed <- lapply(
            names(
                sample_per_enhancer[
                    # filter for enhancers with non-zero tags per sample
                    sample_per_enhancer[,sample_name]>0,][,sample_name]
                ),
            # extract chr and coordinates into separate entities
            function(x){
                unlist(strsplit(x, ":|-"))
                }
        )
        # save list of coordinates to bedfile
        sink(file.path(
            cagerFolder,
            "enhancers",
            paste0(sample_name, "_enhancers.bed")))
        for (line in sample_bed){
            width <- as.integer(line[3]) - as.integer(line[2])
            enhancer_widths[[sample_name]] <- append(
                enhancer_widths[[sample_name]],
                width)
            cat(paste0(line[1], "\t", line[2], "\t", line[3], "\n"))
            }
        sink()
    }
    print("Enhancers per sample bed files saved")
    return(enhancer_widths)
}

plot_enhancer_widths <- function(enhancer_widths, cagerFolder){
  
  enhancer_width_tbl <- tibble(name=character(), width=numeric())

  for (sample in names(enhancer_widths)){
    enhancer_width_tbl <- enhancer_width_tbl %>% add_row(
      name = sample,
      width = enhancer_widths[[sample]])
  }
  
  enhancer_width_plots <- ggplot2::ggplot(data=enhancer_width_tbl, aes(x=width)) +
    geom_histogram(color="grey", fill="grey") +
    facet_wrap(. ~name) +
    ggtitle("Enhancer widths")
  
  tagcluster_width_plots_path <- file.path(cagerFolder, "enhancer_width_plots.pdf")
  ggplot2::ggsave(tagcluster_width_plots_path, enhancer_width_plots)
  return("Sample enhancer widths plotted")
}

save_to_bed <- function(brcage, cagerFolder){
    rtracklayer::export.bedGraph(
        brcage[["enhancers"]]@rowRanges,
        file.path(cagerFolder, "enhancers", "enhancers.bedGraph"))
    return("Enhancers saved to bedGraph file")
}

count_number_of_enhancers <- function(sample_per_enhancer) {
  sample_enhancer_count <- list()
  for (sample in CAGEr::sampleLabels(brcage)) {
    sample_enhancer_count[[sample]] <- sum(as.vector(sample_per_enhancer[,sample])>0)
  }
  return(sample_enhancer_count)
}

plot_pcs <- function(sample_per_enhancer, cagerFolder, pcarank){
  pca_out <- stats::prcomp(sample_per_enhancer, rank. = pcarank)
  plca_to_plot <- as_tibble(data.frame(pca_out$x[,c("PC1","PC2")]), rownames = "name")
  pca_plot <- ggplot2::ggplot(plca_to_plot, aes(x=PC1, y=PC2,label=name)) +
    geom_point() +
    ggrepel::geom_text_repel(hjust=0, vjust=0,size=2,max.overlaps=5)
  ggplot2::ggsave(
    file.path(cagerFolder, "enhancers_pca.pdf"),
    pca_plot
  )
  return("PCA plotted")
}

save_to_file <- function(brcage, sample_per_enhancer, cagerFolder){
    outFileNameSamples <- file.path(
        cagerFolder, "enhancers", "sample_per_enhancer.tsv")
    write.table(sample_per_enhancer, file=outFileNameSamples, quote=FALSE, sep='\t')
    print("Enhancer counts per sample file saved")
    outFileName <- file.path(
        cagerFolder, "enhancers", "enhancers.rds")
    saveRDS(brcage, outFileName)
    return("Enhancers cager rds file saved")
}
