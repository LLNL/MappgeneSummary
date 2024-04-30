## code to prepare `zika_reference` dataset goes here

zika_reference = readxl::read_excel("inst/extdata/zika_KJ776791.xlsx") |>
  dplyr::filter(!SHORT_NAME %in% c("prM", "C"))

zika_genes = zika_reference |>
  ggplot2::ggplot(ggplot2::aes(xmin = GENOME_START,
                               y = "KJ776791",
                               xmax = GENOME_STOP,
                               label = SHORT_NAME,
                               fill = "#037bcf")) +
  gggenes::geom_gene_arrow(arrowhead_height = grid::unit(8, "mm"),
                           arrowhead_width = grid::unit(2, "mm"),
                           arrow_body_height = grid::unit(8, "mm")) +
  ggplot2::geom_text(angle = 90,
            hjust = 0.5, size = 4, ggplot2::aes(x = (GENOME_START + GENOME_STOP) / 2)) +
  gggenes::theme_genes() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_fill_identity() +
  ggplot2::labs(x = "KJ776791 Coordinates")





usethis::use_data(zika_reference, overwrite = TRUE)
usethis::use_data(zika_genes, overwrite = TRUE)

