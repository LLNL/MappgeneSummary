#' Plot trend data by location from outbreak API
#'
#' @param df A tibble or dataframe from one of the get outbreak functions
#' @param title Override the default title
#'
#' @export
#'
#' @return A ggplot2 object which can be modified inline
#'
#' @examples
#' plot_outbreak_trends(top20_results)
#' get_state_outbreak_info(states = c("CA", "OR"), SNVs = "S:Y1155H") %>% plot_outbreak_trends()


plot_outbreak_trends = function(df, title = "Proportions (free y-axis scale)") {

  # guess the number of columns, keeping max of 6 rows
  ncols = ceiling((pull(df, lineage) %>% na.omit() %>% unique() %>% length()) / 6)

  df %>%
    filter(total_count > 0) %>%
    group_by(lineage) %>%
    ggplot(aes(x = date, y = 100 * proportion, fill = location, color = location)) +
      geom_line(aes(group = location), alpha = 0.8) +
      facet_wrap(~lineage, ncol = ncols, scales = "free_y") +
      labs(title = title, y = "Proportion (%)", x = "Month-Year") +
      scale_x_date(date_breaks = "months" , date_labels = "%b-%y")
}
