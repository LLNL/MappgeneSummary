#' Barchart summary of the SNPsift Effects
#'
#' @param df Dataframe or tibble
#'
#' @export
#'

plot_effect_counts <- function(df, keeper_effects = c(
                                 "stop_gained",
                                 "missense_variant",
                                 "conservative_inframe_deletion",
                                 "conservative_inframe_insertion",
                                 "disruptive_inframe_deletion",
                                 "disruptive_inframe_insertion"
                               )) {
  df |>
    dplyr::group_by(EFFECT, PIPELINE) |>
    dplyr::count() |>
    ggplot2::ggplot(ggplot2::aes(x = n, y = EFFECT)) +
    ggplot2::geom_col(
      position = "dodge",
      ggplot2::aes(fill = ifelse(EFFECT %in% keeper_effects, "#6495ed", "gray30"))
    ) +
    ggplot2::geom_text(ggplot2::aes(label = n, x = n), color = "#AAAAAA", hjust = 1) +
    ggplot2::scale_fill_identity() +
    ggplot2::facet_wrap(~PIPELINE) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE), trans = "log10")
}
