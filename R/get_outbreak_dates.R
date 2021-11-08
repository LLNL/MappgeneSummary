#' Get peak dates and proportions from outbreak API data
#'
#' Percents are determined to 0.001% accuracy using the scales package.
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#'
#' @export
#'
#' @examples
#' get_outbreak_peak(top20_results)
#' get_state_outbreak_info(states = c("CA", "OR"), SNVs = "S:Y1155H") %>% get_outbreak_peak()

get_outbreak_peak = function(df) {
  df %>%
    group_by(LOCATION, SNV) %>%
    filter(proportion == max(proportion)) %>%
    mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "PEAK", CODE = paste(LOCATION, SNV, sep = "_")) %>%
    rename("DATE" = "date") %>%
    select(LOCATION, SNV, TYPE, DATE, PERCENT, CODE) %>%
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
#' get_state_outbreak_info(states = c("CA", "OR"), SNVs = "S:Y1155H") %>% get_outbreak_start()

get_outbreak_start= function(df) {

  df %>%
    group_by(LOCATION, SNV) %>%
    filter(lineage_count > 0) %>%
    filter(date == min(date)) %>%
    mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "FIRST", CODE = paste(LOCATION, SNV, sep = "_")) %>%
    rename("DATE" = "date") %>%
    select(LOCATION, SNV, TYPE, DATE, PERCENT, CODE) %>%
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
#' get_state_outbreak_info(states = c("CA", "OR"), SNVs = "S:Y1155H") %>% get_outbreak_latest()

get_outbreak_latest = function(df) {
  df %>%
    group_by(LOCATION, SNV) %>%
    filter(lineage_count > 0) %>%
    filter(date == max(date)) %>%
    mutate(PERCENT = scales::percent(proportion, accuracy = .001), TYPE = "LATEST", CODE = paste(LOCATION, SNV, sep = "_")) %>%
    rename("DATE" = "date") %>%
    select(LOCATION, SNV, TYPE, DATE, PERCENT, CODE) %>%
    unique()
}
