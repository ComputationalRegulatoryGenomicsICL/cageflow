# 
# Parse incoming data from channel
# 

parse_input <- function(sample_info, data_type){
    cleaned_string <- gsub("\\[", "", sample_info[[1]])
    cleaned_string <- gsub("\\]", "", cleaned_string)
    elements <- strsplit(cleaned_string, ",")[[1]]

    # debug
    cat(paste(c("cleaned_string=", cleaned_string), collapse = "|"), "\n")
    cat(paste(c("elements=", elements), collapse = "|"), "\n")
    for(el.i in 1:length(elements)) {
        cat(paste(c("el.i:", el.i, ":'", elements[el.i], "'"),
                  collapse = "|"), "\n")
    }
    # end_debug

    if (tolower(data_type) == "bam") {
        elements <- gsub(" ", "", elements)
        sample_matrix <- matrix(elements, ncol=3, byrow=TRUE)
        sample_table <- as.data.frame(
            sample_matrix, stringsAsFactors=FALSE)
        colnames(sample_table) <- c("id", "single_end", "path")
    } else if (tolower(data_type) == "bigwig") {
        sample_matrix <- matrix(elements, ncol=4, byrow=TRUE)
        sample_table <- as.data.frame(
            sample_matrix, stringsAsFactors=FALSE)
        colnames(sample_table) <- c("id", "single_end", "path1", "path2")
        sample_table$path <- paste(sample_table$path1, sample_table$path2, sep=",")
    } else {
        stop("Either bigwig or bam files should be provided")
    }

    return(sample_table)
}