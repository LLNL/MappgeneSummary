#' Load the mappgene data into a tibble
#'
#' @param path The path to the folder with the .ivar.snpSIFT.txt and .ivar.lofreq.snpSIFT.txt files.
#' @param workflow The mappgene workflow type(s). Accepts a single string, or a vector of both "ivar" and "lofreq".
#' @param sample_cleanup Additional text to remove from the sample names.
#'
#' @export
#'

read_mappgene_data = function(path, workflow = c("ivar", "lofreq"), sample_cleanup = "LLNL_") {

  if (length(setdiff(workflow, c("ivar", "lofreq"))) > 0) {
    stop("Sorry, there are workflow terms I don't recognize. Please use lofreq and/or ivar only")
  }

  l = list()

  if ("ivar" %in% workflow) {
    pattern = "*.ivar.snpSIFT.txt$"

    files = list.files(path = path, pattern = pattern, recursive = TRUE, full.names = T)
    message(paste("Found", length(files), "ivar files"))

    ivar = sapply(files,
                  readr::read_delim,
                  delim = "\t",
                  comment = "##",
                  simplify = FALSE,
                  col_types = "ficcdifffcciic") |>
      lapply(\(x) dplyr::mutate(x, dplyr::across(ALT_QUAL, as.double))) |> # to fix when a file is empty
      dplyr::bind_rows(.id = "id")

    ivar = ivar |>
      dplyr::filter(ALT_QUAL > 20) |>
      dplyr::filter(FILTER == 'PASS')

    l$ivar = ivar
  }

  if ("lofreq" %in% workflow) {
    pattern = "*.ivar.lofreq.snpSIFT.txt$"

    files = list.files(path = path, pattern = pattern, recursive = TRUE, full.names = T)
    message(paste("Found", length(files), "lofreq files"))

    lofreq = sapply(files,
                    readr::read_delim,
                    delim = "\t",
                    comment = "##",
                    simplify = FALSE,
                    col_types = "ficcdifffcciic") |>
      #lapply(\(x) mutate(x, across(ALT_QUAL, as.double))) |> # to fix when a file is empty
      dplyr::bind_rows(.id = "id")
    l$lofreq = lofreq
  }

  df = dplyr::bind_rows(l, .id = "PIPELINE") |>
    dplyr::mutate(FILE = basename(id)) |>
    dplyr::select(-id) |>
    dplyr::mutate(SAMPLE = FILE) |>
    dplyr::mutate(SAMPLE = stringr::str_replace(SAMPLE, ".ivar.snpSIFT.txt", "")) |>
    dplyr::mutate(SAMPLE = stringr::str_replace(SAMPLE, ".ivar.lofreq.snpSIFT.txt", "")) |>
    dplyr::mutate(SAMPLE = stringr::str_replace(SAMPLE, sample_cleanup, ""))

  names(df) = gsub(x = names(df), pattern = "ANN\\[\\*\\]\\.", replacement = "")

  df = df |>
    dplyr::select("SAMPLE", "PIPELINE", "GENE", "FEATUREID", "POS", "REF", "ALT", "AF", "DP", "EFFECT",
           "HGVS_C", "HGVS_P", "CDNA_POS", "AA_POS", "FILE") |>
    unique()|>
    dplyr::mutate(ALT_DP = as.integer(DP * AF)) |>
    tibble::as_tibble()

  return(df)
}
