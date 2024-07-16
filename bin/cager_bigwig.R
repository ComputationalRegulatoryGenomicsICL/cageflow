

args = commandArgs()
bsgenome    = args[6]
bigwig_str  = args[7]
cpus        = args[8]


bigwigs = unlist(stringr::str_split(stringr::str_remove_all(bigwig_str, "[\\[\\],]"), fixed(" ")))

signals = lapply(bigwigs, function(x) rtracklayer::import(x))

names(signals) = basename(bigwigs)

signalsSplit = split(signals, grepl("str1", names(signals)))

plus = lapply(signalsSplit$`TRUE`, function(x) {
  strand(x) = "+"
  return(x)
})

minus = lapply(signalsSplit$`FALSE`, function(x) {
  strand(x) = "-"
  return(x)
})

plus.sample.names = stringr::str_remove_all(names(plus), ".Signal.UniqueMultiple.str1.out.wig.bw")
minus.sample.names = stringr::str_remove_all(names(minus), ".Signal.UniqueMultiple.str2.out.wig.bw")

if (!all(plus.sample.names == minus.sample.names)) {
  if (setequal(plus.sample.names, minus.sample.names)) {
    minus = minus[match(plus.sample.names, minus.sample.names)]
  } else {
    stop("Error: Some basenames of minus- and plus-strand bigWigs are different! Are these bigWigs from different sets of samples? Exit.")
  }
}

merged = mapply(c, plus, minus)

ctss.obj = purrr::reduce(merged, 
                         function(x, y) {
                           full_join(as.data.frame(x), 
                                     as.data.frame(y), 
                                     by = c("seqnames", 
                                            "start", 
                                            "strand")) %>% 
                           dplyr::select(c("seqnames", 
                                           "start", 
                                           "strand"), 
                                         contains("score"))
                         }) %>% 
  dplyr::rename(chr = "seqnames", 
                pos = "start") %>%
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))

names(ctss.obj) = c("chr", "pos", "strand", names(merged))

ctss.obj %<>%
  mutate(across(contains(".bw"), as.integer))

ctss.obj$chr = as.character(ctss.obj$chr)

ctss.obj$strand = as.character(ctss.obj$strand)

ctss.obj = as(ctss.obj, "CAGEexp")

ctss.obj$genomeName = ref.name

rowRanges(ctss.obj@ExperimentList$tagCountMatrix) = as(
  rowRanges(ctss.obj@ExperimentList$tagCountMatrix),
  Class = "CTSS")


