####################
## Plot functions ##
####################

required.libraries <- c(
    "ggplot2",
    "viridis"
)

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

plot_settings <- function(.data, y_value, color_by_value, y_label, title) {
    .data %>% ggplot(aes(
        x = Sample,
        y = {{y_value}},
        colour = {{color_by_value}})) +
    geom_point(size = 5) +
    scale_color_viridis(
        discrete = TRUE,
        begin = 0.1,
        end = 0.9,
        option = "magma") +
    theme_bw(base_size = 8) +
    ylab(y_label) +
    xlab("") + 
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_number_of_ctss <- function(sample_ctss_count, yaxistitle, mytitle, myfilename) {
  sample_ctss_count_table <- tibble::enframe(
    sample_ctss_count,
    name = "Sample",
    value = "count") %>%
    tidyr::unnest(cols=count)
  ctss_plot <- sample_ctss_count_table %>%
    plot_settings(
      y_value = count,
      color_by_value = Sample,
      y_label = yaxistitle,
      title = mytitle)
  return(ctss_plot)
}

plot_pcs <- function(ce_tmp, pcarank){
  pca_out <- prcomp(
    ce_tmp,
    rank. = pcarank)
  plca_to_plot <- as_tibble(
    data.frame(
      pca_out$rotation[,c("PC1","PC2")]),
      rownames = "name")
  pca_plot <- ggplot2::ggplot(
    plca_to_plot,
    aes(x=PC1, y=PC2,label=name)) +
    geom_point() +
    ggrepel::geom_text_repel(
      hjust=0,
      vjust=0,
      size=2,
      max.overlaps=5)

  return(pca_plot)
}
