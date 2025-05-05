save_plot <- function(filename, plot_out){
    filename = paste0("plots/", filename)
    ggsave(
        filename = filename,
        plot = plot_out,
        width = 20,
        height = 10,
        limitsize = FALSE)
    datafilename <- gsub("pdf", "rds", filename)
    saveRDS(plot_out, datafilename)
}