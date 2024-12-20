#' Get peak dates and proportions from outbreak API data
#'
#' If multiple dates have the same proportion, only the fist date is shown.
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export

get_outbreak_peak <- function(df) {
  df |>
    dplyr::group_by(location, lineage) |>
    dplyr::filter(proportion == max(proportion)) |>
    dplyr::filter(date == min(date)) |> # get first date if multiple are the same
    dplyr::mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "PEAK", CODE = paste(location, lineage, sep = "_")) |>
    dplyr::rename("DATE" = "date") |>
    dplyr::select(location, lineage, TYPE, DATE, PERCENT, CODE) |>
    unique()
}

#' Get start dates and proportions from outbreak API data
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export

get_outbreak_start <- function(df) {
  df |>
    dplyr::group_by(location, lineage) |>
    dplyr::filter(lineage_count_rolling > 0) |>
    dplyr::filter(date == min(date)) |>
    dplyr::mutate(
      PERCENT = scales::percent(proportion, accuracy = .001),
      TYPE = "FIRST",
      CODE = paste(location, lineage, sep = "_")
    ) |>
    dplyr::rename("DATE" = "date") |>
    dplyr::select(location, lineage, TYPE, DATE, PERCENT, CODE) |>
    unique()
}

#' Get most recent dates and proportions from outbreak API data
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export

get_outbreak_latest <- function(df) {
  df |>
    dplyr::mutate(date = as.Date(date)) |>
    dplyr::group_by(location, lineage) |>
    dplyr::filter(lineage_count_rolling > 0) |>
    dplyr::filter(date == max(date)) |>
    dplyr::mutate(
      PERCENT = scales::percent(proportion, accuracy = 0.001),
      TYPE = "LATEST",
      CODE = paste(location, lineage, sep = "_")
    ) |>
    dplyr::rename(DATE = "date") |>
    dplyr::select(location, lineage, TYPE, DATE, PERCENT, CODE) |>
    unique()
}
