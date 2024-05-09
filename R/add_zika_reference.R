#' Add data for processed Zika peptides
#'
#' Given a mappgene dataframe, this function will add the genome reference data
#' from Zika KJ776791.
#'
#' @param df A mappgene summary dataframe
#' @param remove_seqs Remove the nucleotide and amino acid sequences from the dataframe
#'
#' @export
#' @importFrom rlang .data


add_zika_reference <- function(df,
                               remove_seqs = TRUE) {
  df <- dplyr::mutate(df, dummy_join = "ZIKA")
  zika_KJ776791 <- dplyr::mutate(MappgeneSummary::zika_reference, dummy_join = "ZIKA")

  # remove pre proteins

  zika_KJ776791 <- zika_KJ776791 |>
    dplyr::filter(!.data$SHORT_NAME %in% c("prM", "C"))

  df <- df |>
    dplyr::inner_join(zika_KJ776791,
      by = "dummy_join",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(.data$POS >= .data$GENOME_START, .data$POS <= .data$GENOME_STOP) |>
    dplyr::select(-c(.data$GENOME_START, .data$GENOME_STOP, .data$dummy_join)) |>
    dplyr::mutate(AA_PROCESSED_POS = .data$AA_POS - .data$AA_START + 1)

  if (remove_seqs) {
    df <- df |>
      dplyr::select(-.data$NT, -.data$AA)
  }

  df
}
