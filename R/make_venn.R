#' Make Pipeline Venn
#'
#' @param df Dataframe or tibble
#' @param slices Column name to fill the slices of the venn
#'
#' @export
#'

make_venn = function(df, slices = "SNV") {

  require(ggVennDiagram)

  lofreq = as_tibble(df) %>%
    filter(PIPELINE == 'lofreq')

  ivar = as_tibble(df) %>%
    filter(PIPELINE== 'ivar')

  x1 = list(lofreq = unique(lofreq[[slices]]),
            ivar = unique(ivar[[slices]]))

  p = ggVennDiagram(x1) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    scale_color_manual(values = c("black", "black")) +
    theme(legend.position = "none") +
    labs(title = slices)

  l = list("data" = x1, "plot" = p)
}
