#' Make coverage files from folder with bedgraphs
#'
#' @param file_path A directory with bedgraph files
#'
#' @export
#' @importFrom rlang .data


get_coverages <- function(file_path) {
  bed_list <- list.files(path = file_path, pattern = "*ivar.bedgraph", full.names = T, recursive = TRUE)

  coverage <- sapply(bed_list, readr::read_delim,
                     col_names = c("C", "START", "STOP", "COV"),
                     delim = "\t",
                     simplify = FALSE,
                     col_types = readr::cols()) |>
    dplyr::bind_rows(.id = "SAMPLE") |>
    dplyr::mutate(FILE = dirname(.data$SAMPLE)) |>
    dplyr::mutate(FILE = stringr::str_remove(.data$FILE, file_path)) |>
    dplyr::mutate(SAMPLE = basename(.data$SAMPLE)) |>
    dplyr::mutate(SAMPLE = stringr::str_remove(.data$SAMPLE, ".ivar.bedgraph")) |>
    dplyr::mutate(POS = purrr::map2(.data$START, .data$STOP - 1, ~ seq(from = .x, to = .y))) |>
    dplyr::mutate(FOLDER = dirname(.data$FILE)) |>
    dplyr::select(.data$SAMPLE, .data$COV, .data$POS, .data$FOLDER) |>
    tidyr::unnest(cols = c(.data$POS)) |>
    dplyr::arrange(.data$POS)

  coverage
}



#' Plot ECDF coverage
#'
#' @param coverage A coverage from dataframe
#' @param pseudocount A pseudocount to add to coverage to make log10 transformation work
#' @param intercept An x value to draw a vertical line
#'
#' @export
#' @importFrom rlang .data

plot_coverage_ecdf <- function(coverage, pseudocount = 1, intercept = 100) {
  coverage |>
    dplyr::mutate(COV = .data$COV + pseudocount) |>
    ggplot2::ggplot(ggplot2::aes(x = COV)) +
    ggplot2::stat_ecdf(pad = F) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = intercept, linetype = 2, color = "#ff0000aa", linewidth = 1) +
    ggplot2::facet_wrap(~SAMPLE) +
    ggplot2::scale_x_continuous(
      labels = function(x) {
        format(x,
          big.mark = ",",
          decimal.mark = ".",
          scientific = FALSE
        )
      },
      trans = "log10"
    ) +
    ggplot2::labs(
      x = "Coverage (log10)",
      y = "Cumulative Distribution"
    )
}
