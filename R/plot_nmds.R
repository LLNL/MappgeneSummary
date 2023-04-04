#' Plot NMDS of AF values
#'
#' If you don't get convergent, try changing the seed with set.seed()
#'
#' @param df Dataframe or tibble
#' @param points The main values to ordinate
#' @param species The secondary values to order the points with
#' @param distance A distance method
#'
#' @export
#'

plot_nmds = function(df, points = "SAMPLE", species = "SNV", distance = "bray") {
  require(vegan)
  require(ggrepel)

  if (distance == "jaccard") {
    set_binary = T
  } else {
    set_binary = F
  }

  df.v1 = df %>%
    select(points, species, AF) %>%
    group_by_at(c(points, species)) %>%
    summarize(AF = mean(AF), .groups = 'drop') %>%
    pivot_wider(names_from = points, values_from = AF, values_fill = list(AF = 0)) %>%
    column_to_rownames(species) %>%
    as.matrix() %>%
    t() %>%
    vegan::metaMDS(distance = distance, binary = set_binary, trace = 0)

  points.df = vegan::scores(df.v1)$sites %>%
    as.data.frame() %>%
    rownames_to_column(points)

  species.df = vegan::scores(df.v1)$species %>%
    as.data.frame() %>%
    rownames_to_column(species)

  g = ggplot(points.df, aes(x = NMDS1, y = NMDS2)) +
    ggrepel::geom_text_repel(aes_string(label = points), size = 3) +
    geom_point(data = species.df, size = 1, alpha = 0.3, aes(x = NMDS1, y= NMDS2)) +
    geom_point(size = 3, pch = 24, alpha = 0.8, fill = "#037bcf") +
    theme(legend.position = "right") +
    labs(title = paste(points, "*", species), subtitle = paste("Distance =", distance, ", Stress =", df.v1$stress * 100))

  return(list("plot" = g, "nmds" = df.v1))
}
