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


#' Get the first, peak and most recent times and plots for a list of mutations
#'
#' @param df A tibble or dataframe from get_outbreak_info
#' @param snv_list A list of the SNVs
#' @param ncol The number of columns to facet the plots
#' @param min.date Filter the df to only include dates later than or equal to this date (YYYY-MM-DD)
#' @param max.date Filter the df to only include dates earlier than or equal to this date (YYYY-MM-DD)
#'
#' @export
#'

plot_outbreak_milestones = function(df, snv_list, ncol = 2, min.date = NA, max.date = NA) {

  require(lubridate)

  if (!is.na(max.date)) {
    max.date = lubridate::ymd(max.date)
    df = df %>%
      filter(date <= max.date)
  }

  if (!is.na(min.date)) {
    min.date = lubridate::ymd(min.date)
    df = df %>%
      filter(date >= min.date)
  }

  snv_list = as.factor(snv_list)

  p = df %>%
    get_outbreak_peak() %>%
    mutate(lineage = as.factor(lineage)) %>%
    mutate(TYPE = as.factor(TYPE)) %>%
    mutate(location = as.factor(location))

  f = df %>%
    get_outbreak_start() %>%
    mutate(lineage = as.factor(lineage)) %>%
    mutate(TYPE = as.factor(TYPE)) %>%
    mutate(location = as.factor(location))

  l = df %>%
    get_outbreak_latest() %>%
    mutate(lineage = as.factor(lineage)) %>%
    mutate(TYPE = as.factor(TYPE)) %>%
    mutate(location = as.factor(location))

  b1 = df %>%
    mutate(lineage = fct_relevel(lineage, levels(snv_list))) %>%
    plot_outbreak_trends() +
      scale_color_manual(values = palette_jak$bay(3)) +
      facet_wrap(~lineage, ncol = ncol, scales = "free_y")

  b2 = bind_rows(f, l, p) %>%
    mutate(lineage = fct_relevel(lineage, levels(snv_list))) %>%
    ggplot(aes(x = DATE, y = location, fill = TYPE)) +
      geom_point(pch = 21, size = 3, alpha = 0.7) +
      scale_fill_manual(values = palette_jak$bay(3)) +
      facet_wrap(~lineage, ncol = ncol) +
      labs(title = "First, last and peak occurrence", x = "Month-Year")

  b3 = bind_rows(f, l, p) %>%
    mutate(TYPE = fct_relevel(TYPE, 'FIRST', 'PEAK', 'LATEST')) %>%
    arrange(TYPE) %>% # some wonkiness to get the line to plot in the right order when two TYPEs have the same data
    mutate(lineage = fct_relevel(lineage, levels(snv_list))) %>%
    mutate(location = fct_relevel(location, 'Worldwide', 'United States', 'California')) %>%
    ggplot(aes(x = DATE, y = TYPE, color = location)) +
      geom_path(aes(group = location), size = 1, alpha = 0.5) + # some wonkiness to get the line to plot in the right order when two TYPEs have the same data
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = palette_jak$bay(3)) +
      facet_wrap(~lineage, ncol = ncol) +
      labs(title = "First, peak and last occurrence", x = "Month-Year")

  l = list("A" = b1, "B" = b2, "C" = b3, "peak" = p, "start" = f, "latest" = l)

}
