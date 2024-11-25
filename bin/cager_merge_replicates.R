# merge replicates - automatic extraction of replicate counts
# only if they follow sampleName_repX pattern where X is ordered and starts at the same number for all samples (eg 1)

merge_replicates <- function(ce){
    if (!all(grepl( "_rep", sampleLabels(ce), fixed = TRUE))){
        print("No replicates with syntax _repX is found in all samples, no merging is performed")
        return(ce)
    }
    # extract sample names
    sample_names <- unlist(unique(lapply(
        sampleLabels(ce),
        function(x) {
            unlist(strsplit(x, "_rep"))[1]
        }
    )))
    # extract replicate number per sample
    sample_replicate_counts <- unlist(lapply(
        sampleLabels(ce),
        function(x) {
            unlist(strsplit(x, "_rep"))[2]
        }
    ))
    # count how many replicates there are of each sample
    # assumes that each sample is names <samplename>_rep<replicatenumber> 
    # and they are ordered when inserted (ordering is done before reading in data to CAGEr)
    new_indeces <- c()
    index_count <- 0
    index_location <- 1
    # it can happen that the replicate counts do not start at 1, but at 3 or something
    # if there are replicates 1, 2, then 3 and 4 for a different sample, this will not work properly
    minimum_rep_count = min(sample_replicate_counts)
    for (src in sample_replicate_counts){
        if (src == minimum_rep_count){
            index_count <- index_count + 1
        }
        new_indeces[index_location] <- index_count
        index_location <- index_location + 1
    }
    # merge replicates
    ce <- mergeSamples(
        ce,
        mergeIndex = new_indeces,
        mergedSampleLabels = sample_names
    )
    return(ce)
}