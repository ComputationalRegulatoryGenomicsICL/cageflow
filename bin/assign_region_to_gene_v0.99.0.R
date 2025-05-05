#' Assign Region to the Promoter of the Closest Gene
#'
#' This function assigns a region to the promoter of the closest gene using a 
#' specified upstream and downstream extension from the transcription start 
#' site (TSS).
#'
#' @param region A `GRanges` object representing the input regions.
#' @param txdb A TxDb object containing the gene annotations.
#' @param upstream An integer specifying the number of base pairs upstream of 
#' the TSS to include in the promoter region (default: 500).
#' @param downstream An integer specifying the number of base pairs downstream 
#' of the TSS to include in the promoter region. The downstream extension will 
#' be limited to the end of the transcript (default: 500).
#' @return A `GRanges` object with the original regions and the closest gene 
#' information, including gene ID, gene coordinates, strand, and distance to 
#' the promoter.
#' @import GenomicRanges
#' @import AnnotationDbi
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene transcripts promoters
#' @export
#' @examples
#' # Example usage:
#' # Define a sample region as a GRanges object
#' sample_region <- GRanges(seqnames = "chr1", 
#'                          ranges = IRanges(start = 1000000, end = 1000500))
#'
#' # Load TxDb object (example with human genome)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' # Call the function with the sample region
#' result <- assign_region_to_promoter(sample_region, 
#'                                     TxDb.Hsapiens.UCSC.hg38.knownGene)
#' print(result)

assign_region_to_promoter <- function(region, txdb, upstream = 500,
                                      downstream = 500) {
  # Extract transcripts from the TxDb object
  transcripts <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
  
  # Extend promoters using the promoters function
  # Adjust the downstream end to the minimum of the downstream value and the 
  # end of the transcript

  promoters_extended <- promoters(transcripts, upstream = upstream, 
                                  downstream = pmin(downstream, 
                                                    width(transcripts)))
  
  # end(promoters_extended) <- pmin(end(promoters_extended), 
  #                                 end(transcripts[match(
  #                                       names(promoters_extended), 
  #                                       names(transcripts))]))
  
  # Find overlaps between the input regions and the extended promoters
  overlaps <- findOverlaps(region, promoters_extended)
  
  # Get the closest promoter for each region
  closest_hits <- as.data.frame(overlaps)
  # closest_hits$distance <- abs(
  #   start(region[closest_hits$queryHits]) - 
  #     start(promoters_extended[closest_hits$subjectHits]))
  closest_hits$distance <- distance(region[queryHits(overlaps)], 
                                    resize(transcripts[subjectHits(overlaps)],
                                           width = 1))
  
  # If a region overlaps multiple promoters, take the closest one
  closest_hits <- closest_hits[order(closest_hits$queryHits, 
                                     closest_hits$distance), ]
  closest_hits <- closest_hits[!duplicated(closest_hits$queryHits), ]
  
  # Retrieve the corresponding gene IDs for the closest promoters
  closest_genes <- transcripts[closest_hits$subjectHits]
  
  # Create a data frame with the original regions and the closest gene information
  result <- data.frame(
    region_seqnames = seqnames(region)[closest_hits$queryHits],
    region_start = start(region)[closest_hits$queryHits],
    region_end = end(region)[closest_hits$queryHits],
    mcols(region[closest_hits$queryHits]),
    closest_gene_id = unlist(closest_genes$gene_id),
    closest_gene_seqnames = seqnames(closest_genes),
    closest_gene_start = start(closest_genes),
    closest_gene_end = end(closest_genes),
    closest_gene_strand = strand(closest_genes),
    distance_to_promoter = closest_hits$distance
  ) |> makeGRangesFromDataFrame(seqnames.field = "region_seqnames",
                                  start.field = "region_start",
                                  end.field = "region_end",
                                  keep.extra.columns = TRUE)
  
  return(result)
}
