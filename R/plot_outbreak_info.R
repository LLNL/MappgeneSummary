#' Plot trend data by location from outbreak API
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#' @param title Override the default title
#'
#' @export
#'
#' @return A ggplot2 object which can be modified inline

plot_outbreak_trends <- function(df, title = "Proportions (free y-axis scale)", line_width = 1, line_alpha = 0.8) {
  # guess the number of columns, keeping max of 6 rows
  ncols <- ceiling((dplyr::pull(df, lineage) |> na.omit() |> unique() |> length()) / 6)

  df |>
    dplyr::filter(total_count > 0) |>
    dplyr::group_by(lineage) |>
    ggplot2::ggplot(ggplot2::aes(x = date, y = 100 * proportion, fill = location, color = location)) +
    ggplot2::geom_line(ggplot2::aes(group = location), alpha = line_alpha, size = line_width) +
    ggplot2::facet_wrap(~lineage, ncol = ncols, scales = "free_y") +
    ggplot2::labs(title = title, y = "Proportion (%)", x = "Month-Year") +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y")
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

plot_outbreak_milestones <- function(df, snv_list, ncol = 2, min.date = NA, max.date = NA, colors = c(), line_width = 1, line_alpha = 0.8) {
  require(lubridate)

  if (!is.na(max.date)) {
    max.date <- lubridate::ymd(max.date)
    df <- df |>
      dplyr::filter(date <= max.date)
  }

  if (!is.na(min.date)) {
    min.date <- lubridate::ymd(min.date)
    df <- df |>
      dplyr::filter(date >= min.date)
  }

  snv_list <- as.factor(snv_list)

  # color palette
  locations <- as.factor(unique(df$location))

  color_palette <- grDevices::colorRampPalette(c("#00496f", "#0f85a0", "#edd746", "#ed8b00", "#dd4124"))

  if (length(colors) == 0) {
    colors <- color_palette(length(locations))
  }

  p <- df |>
    get_outbreak_peak() |>
    dplyr::mutate(lineage = as.factor(lineage)) |>
    dplyr::mutate(TYPE = as.factor(TYPE)) |>
    dplyr::mutate(location = as.factor(location))

  f <- df |>
    get_outbreak_start() |>
    dplyr::mutate(lineage = as.factor(lineage)) |>
    dplyr::mutate(TYPE = as.factor(TYPE)) |>
    dplyr::mutate(location = as.factor(location))

  l <- df |>
    get_outbreak_latest() |>
    dplyr::mutate(lineage = as.factor(lineage)) |>
    dplyr::mutate(TYPE = as.factor(TYPE)) |>
    dplyr::mutate(location = as.factor(location))

  b1 <- df |>
    dplyr::mutate(lineage = fct_relevel(lineage, levels(snv_list))) |>
    plot_outbreak_trends(line_width = line_width, line_alpha = line_alpha) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::facet_wrap(~lineage, ncol = ncol, scales = "free_y") +
    ggplot2::theme(legend.position = "right")

  b2 <- dplyr::bind_rows(f, l, p) |>
    dplyr::mutate(lineage = fct_relevel(lineage, levels(snv_list))) |>
    ggplot2::ggplot(ggplot2::aes(x = DATE, y = location, shape = TYPE)) +
    ggplot2::geom_point(size = 3) +
    # scale_fill_manual(values = color_palette(3)) +
    ggplot2::scale_shape_manual(values = c("FIRST" = 3, "PEAK" = 2, "LATEST" = 4)) +
    ggplot2::facet_wrap(~lineage, ncol = ncol) +
    ggplot2::labs(title = "First, last and peak occurrence", x = "Month-Year") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y")

  l <- list("A" = b1, "B" = b2, "peak" = p, "start" = f, "latest" = l)
}
