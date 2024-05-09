#' Plot Heatmap of Specific Gene
#'
#' If working with a lot of data it is recommended to filter your df first, and then pass it into the plot_heatmap() function
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#' @param columns The column name from df to use as columns in the heatmap
#' @param rows The column name from df to use as rows in the heatmap
#' @param na.value The numerical value to set the NA values to. Changing this (from 0 to -1) can affect the "chaining" effect of clustering.
#'
#' @export
#'
#' @examples
#' plot_heatmap(df, rows = "SNV", columns = "SAMPLE")
#' df |> filter(GENE == "ORF3a") |> plot_heatmap(rows = "EFFECT", columns = "SAMPLE")

plot_heatmap = function(df, columns = "SAMPLE", rows = "SNV", na.value = 0) {

  df |>
    tibble::as_tibble() |>
    dplyr::select(columns, rows, AF) |>
    dplyr::group_by_at(c(columns, rows)) |>
    dplyr::summarize(AF = mean(AF), .groups = 'drop') |>
    tidyr::pivot_wider(names_from = columns, values_from = AF, values_fill = (AF = na.value)) |>
    tibble::column_to_rownames(rows) |>
    as.matrix() |>
    pheatmap::pheatmap(show_rownames = T,
                       border_color = NA,
                       color = viridis::viridis(100))
}
