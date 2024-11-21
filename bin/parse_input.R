# 
# Parse incoming data from channel
# 

parse_input <- function(sample_info){
    sample_table = read.delim(
        sample_info,
        header = TRUE,
        sep = ",")

    return(sample_table)
}