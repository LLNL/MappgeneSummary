## code to prepare `zika_reference` dataset goes here

zika_reference = readxl::read_excel("inst/extdata/zika_KJ776791.xlsx") |>
  dplyr::filter(!SHORT_NAME %in% c("prM", "C"))

usethis::use_data(zika_reference, overwrite = TRUE)
