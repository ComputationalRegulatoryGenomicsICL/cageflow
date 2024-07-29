# helper functions copied from https://github.com/charles-plessy/CAGEr/blob/master/R/Annotations.R
ranges2names <- function(rangesA, rangesB) {
  if (is.null(names(rangesB)))
    stop(sQuote("rangesB"), " must have names.")
  names <- findOverlaps(rangesA, rangesB)
  names <- as(names, "List")
  names <- extractList(names(rangesB), names)
  names <- unique(names)
  names <- unstrsplit(names, ";")
  Rle(names)
}

ranges2genes <- function(ranges, genes) {
  if (is.null(genes$gene_name))
    stop("Annotation must contain ", dQuote("gene_name"), " metdata.")
  names(genes) <- genes$gene_name
  ranges2names(ranges, genes)
}

# helper function copied from https://github.com/charles-plessy/CAGEr/commit/23473dfbad45cb645d3143f12ba47354764318d6
txdb2annot <- function(ranges, annot) {
  findOverlapsBool <- function(A, B) {
    overlap <- findOverlaps(A, B)
    overlap <- as(overlap, "List")
    any(overlap)
  }

  classes <- c("promoter", "exon", "intron", "unknown")
  p <- findOverlapsBool(ranges, trim(suppressWarnings(promoters(annot, 500, 500))))
  e <- findOverlapsBool(ranges, exons(annot))
  t <- findOverlapsBool(ranges, transcripts(annot))
  annot <- sapply( 1:length(ranges), function(i) {
    if      (p[i]) {classes[1]}
    else if (e[i]) {classes[2]}
    else if (t[i]) {classes[3]}
    else           {classes[4]}
  })

  annot <- factor(annot, levels = classes)
  Rle(annot)
}

# additional helper functions
annotate_gene_regions <- function(ce, annot_data, debugMode){
    annot_genes <- genes(annot_data)
    annot_genes$gene_name <- annot_genes$gene_id
    annotateCTSS(ce, annot_genes)
    CTSScoordinatesGR(ce)$genes      <- ranges2genes(CTSScoordinatesGR(ce), annot_genes)
    CTSScoordinatesGR(ce)$annotation <- txdb2annot(CTSScoordinatesGR(ce), annot_data)
    annot <- sapply(
        CTSStagCountDF(ce),
        function(X) tapply(
            X,
            CTSScoordinatesGR(ce)$annotation,
            sum
            )
        )
    colData(ce)[levels(CTSScoordinatesGR(ce)$annotation)] <- DataFrame(t(annot))

    if (debugMode){
        print(validObject(ce))
        print(colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")])
    }
    return(ce)
}
