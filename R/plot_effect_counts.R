#' Barchart summary of the SNPsift Effects
#'
#' @param df Dataframe or tibble
#'
#' @export
#'

plot_effect_counts = function(df, keeper_effects = c('stop_gained',
                                                     'missense_variant',
                                                     'conservative_inframe_deletion',
                                                     'conservative_inframe_insertion',
                                                     'disruptive_inframe_deletion',
                                                     'disruptive_inframe_insertion')) {

  df %>%
    group_by(EFFECT, PIPELINE) %>%
    count() %>%
    ggplot(aes(x = n, y = EFFECT)) +
      geom_col(position = "dodge", aes(fill = ifelse(EFFECT %in% keeper_effects, "#6495ed", "gray30"))) +
      geom_text(aes(label = n, x = n), color = "#AAAAAA", hjust = 1) +
      scale_fill_identity() +
      facet_wrap(~PIPELINE) +
      theme(legend.position = "none") +
      scale_x_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE), trans = "log10")
}
