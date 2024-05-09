#' Add a column with SNV data
#'
#' @param df A dataframe of imported snpSIFT data from mappgene output
#'
#' @export
#' @importFrom rlang .data

add_SNV <- function(df,
                    gene_column = "SHORT_NAME") {

  if (!gene_column %in% colnames(df)) {
    warning(glue::glue("Column {gene_column} not found in dataframe. Using default column name GENE"))
    gene_column <- "GENE"
  }

  df |>
    dplyr::arrange(.data$POS) |>
    dplyr::mutate(SNV = paste(!!as.name(gene_column), .data$POS, .data$REF, .data$ALT, sep = "."))
}
