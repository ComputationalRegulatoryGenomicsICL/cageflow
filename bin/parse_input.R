# 
# Parse incoming data from channel
# 

parse_input <- function(sample_info){
    sample_table = read.delim(
        sample_info,
        header = FALSE,
        sep = "\t")
    names(sample_table) = c("id", "single_end", "path")

    return(sample_table)
}