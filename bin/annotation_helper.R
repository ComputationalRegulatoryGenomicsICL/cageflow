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
annotate_gene_regions <- function(brcage, annot_data, debugMode){
    annot_genes <- genes(annot_data)
    annot_genes$gene_name <- annot_genes$gene_id
    annotateCTSS(brcage, annot_genes)
    CTSScoordinatesGR(brcage)$genes      <- ranges2genes(CTSScoordinatesGR(brcage), annot_genes)
    CTSScoordinatesGR(brcage)$annotation <- txdb2annot(CTSScoordinatesGR(brcage), annot_data)
    annot <- sapply(
        CTSStagCountDF(brcage),
        function(X) tapply(
            X,
            CTSScoordinatesGR(brcage)$annotation,
            sum
            )
        )
    colData(brcage)[levels(CTSScoordinatesGR(brcage)$annotation)] <- DataFrame(t(annot))

    if (debugMode){
        print(validObject(brcage))
        print(colData(brcage)[,c("librarySizes", "promoter", "exon", "intron", "unknown")])
    }
    return(brcage)
}

# workaround so that txdb annotation works with consensus
setMethod("annotateConsensusClusters", c("CAGEexp", "TxDb"), function (object, ranges){
  if(is.null(experiments(object)$tagCountMatrix))
    stop("Input does not contain CTSS expressiond data, see ", dQuote("getCTSS()"), ".")
  consensusClustersGR(object)$annotation <- txdb2annot(consensusClustersGR(object), ranges)
  consensusClustersGR(object)$genes    <- ranges2genes(consensusClustersGR(object), ranges)
  validObject(object)
  object
})

setMethod("annotateCTSS", c("CAGEexp", "TxDb"), function (object, ranges){
  g <- genes(ranges)
  g$gene_name <- g$gene_id
  annotateCTSS(object, g)
  CTSScoordinatesGR(object)$genes      <- ranges2genes(CTSScoordinatesGR(object), g)
  CTSScoordinatesGR(object)$annotation <- txdb2annot(CTSScoordinatesGR(object), ranges)

  annot <- sapply( CTSStagCountDF(object)
                   , function(X) tapply(X, CTSScoordinatesGR(object)$annotation, sum))
  colData(object)[levels(CTSScoordinatesGR(object)$annotation)] <- DataFrame(t(annot))

  validObject(object)
  object
})