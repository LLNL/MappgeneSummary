#' Make Pipeline Venn
#'
#' @param df Dataframe or tibble
#' @param slices Column name to fill the slices of the venn
#'
#' @export
#'

make_venn <- function(df, slices = "SNV") {

  lofreq <- tibble::as_tibble(df) |>
    dplyr::filter(PIPELINE == "lofreq")

  ivar <- tibble::as_tibble(df) |>
    dplyr::filter(PIPELINE == "ivar")

  x1 <- list(
    lofreq = unique(lofreq[[slices]]),
    ivar = unique(ivar[[slices]])
  )

  p <- ggVennDiagram::ggVennDiagram(x1) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::scale_color_manual(values = c("black", "black")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = slices)

  l <- list("data" = x1, "plot" = p)
}
