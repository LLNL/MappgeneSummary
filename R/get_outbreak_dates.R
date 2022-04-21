#' Get peak dates and proportions from outbreak API data
#'
#' If multiple dates have the same proportion, only the fist date is shown.
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export
#'
#' @examples
#' get_outbreak_peak(top20_results)
#' get_outbreak_info(SNVs = "S:Y1155H") %>% get_outbreak_peak()

get_outbreak_peak = function(df) {

  df %>%
    group_by(location, lineage) %>%
    filter(proportion == max(proportion)) %>%
    filter(date == min(date)) %>% # get first date if multiple are the same
    mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "PEAK", CODE = paste(location, lineage, sep = "_")) %>%
    dplyr::rename("DATE" = "date") %>%
    select(location, lineage, TYPE, DATE, PERCENT, CODE) %>%
    unique()
}

#' Get start dates and proportions from outbreak API data
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export
#'
#' @examples
#' get_outbreak_start(top20_results)
#' get_outbreak_info(SNVs = "S:Y1155H") %>% get_outbreak_start()

get_outbreak_start= function(df) {

  df %>%
    group_by(location, lineage) %>%
    filter(lineage_count_rolling > 0) %>%
    filter(date == min(date)) %>%
    mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "FIRST", CODE = paste(location, lineage, sep = "_")) %>%
    dplyr::rename("DATE" = "date") %>%
    select(location, lineage, TYPE, DATE, PERCENT, CODE) %>%
    unique()
}

#' Get most recent dates and proportions from outbreak API data
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export
#'
#' @examples
#' get_outbreak_latest(top20_results)
#' get_outbreak_info(SNVs = "S:Y1155H") %>% get_outbreak_latest()

get_outbreak_latest = function (df) {
  df %>%
    mutate(date = as.Date(date)) %>%
    group_by(location, lineage) %>%
    filter(lineage_count_rolling > 0) %>%
    filter(date == max(date)) %>%
    mutate(PERCENT = scales::percent(proportion,accuracy = 0.001), TYPE = "LATEST", CODE = paste(location, lineage, sep = "_")) %>%
    dplyr::rename(DATE = "date") %>%
    select(location, lineage, TYPE, DATE, PERCENT, CODE) %>%
    unique()
}



