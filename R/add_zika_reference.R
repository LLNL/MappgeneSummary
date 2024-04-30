#' Add data for processed Zika peptides
#'
#' Given a mappgene dataframe, this function will add the genome reference data
#' from Zika KJ776791.
#'
#' @param df A mappgene summary dataframe
#' @param remove_seqs Remove the nucleotide and amino acid sequences from the dataframe
#'
#' @export


add_zika_reference <- function(df,
                               remove_seqs = TRUE) {
  df <- dplyr::mutate(df, dummy_join = "ZIKA")
  zika_KJ776791 <- dplyr::mutate(zika_reference, dummy_join = "ZIKA")

  # remove pre proteins

  zika_KJ776791 <- zika_KJ776791 |>
    dplyr::filter(!SHORT_NAME %in% c("prM", "C"))

  df <- df |>
    dplyr::inner_join(zika_KJ776791,
      by = "dummy_join",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(POS >= GENOME_START, POS <= GENOME_STOP) |>
    dplyr::select(-c(GENOME_START, GENOME_STOP, dummy_join)) |>
    dplyr::mutate(AA_PROCESSED_POS = AA_POS - AA_START + 1)

  if (remove_seqs) {
    df <- df |>
      dplyr::select(-NT, -AA)
  }

  df
}
