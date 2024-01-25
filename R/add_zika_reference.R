#' Add data for processed Zika peptides
#'
#' @param df A mappgene summary dataframe
#'
#' @export


add_zika_reference = function(df) {

  df = mutate(df, dummy_join = "ZIKA")
  zika_KJ776791 = mutate(zika_reference, dummy_join = "ZIKA")

  # remove pre proteins

  zika_KJ776791 = zika_KJ776791 |>
    filter(!SHORT_NAME %in% c("prM", "C"))

  df |>
    inner_join(zika_KJ776791,
               by = "dummy_join",
               relationship = "many-to-many") |>
    filter(POS >= GENOME_START, POS <= GENOME_STOP ) |>
    select(-c(GENOME_START, GENOME_STOP, dummy_join)) |>
    mutate(AA_PROCESSED_POS = AA_POS - AA_START + 1)

}
