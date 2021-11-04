#' Copy Summary Rmd Template
#'
#' @param path The path for the .Rmd and style.css files
#' @param name The name for the .Rmd file
#'
#' @export

copy_summary_template = function(path, name = "summary.Rmd", overwrite = TRUE) {
  file.copy(system.file("Rmd", "summary_template.Rmd", package = "MappgeneSummary"), file.path(path, name), overwrite = overwrite)
  file.copy(system.file("Rmd", "style.css", package = "MappgeneSummary"), path, overwrite = overwrite)
}
