#' Parse REF and ALT AA Characters
#'
#' The HGVS_P is parsed according to 11 different rules, depending on the EFFECT,
#' and if there is a deletion, insertion or duplication.
#'
#' @param df A dataframe of imported snpSIFT data from mappgene output
#'
#' @export

parse_ref_alt = function(df) {

  patterns = c("p\\.", "Phe", "Leu", "Ile", "Met", "Val", "Ser", "Pro", "Thr",
               "Ala", "Tyr", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys",
               "Trp", "Arg", "Gly")

  replacement = c("", "F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "H",
  "Q", "N", "K", "D", "E", "C", "W", "R", "G")

  df = df %>%
    mutate(HGVS_P_modified = HGVS_P)

  a = df %>%
    filter(EFFECT == 'missense_variant') %>%
    mutate(HGVS_P_modified = str_replace_all(HGVS_P_modified, setNames(replacement, patterns))) %>%
    separate(HGVS_P_modified, sep = "[0-9]+", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    filter(REF_AA != ALT_AA)

  b = df %>%
    filter(EFFECT %in% c('conservative_inframe_deletion', 'disruptive_inframe_deletion')) %>%
    separate(HGVS_P_modified, sep = "del", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    mutate(REF_AA = str_replace_all(REF_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = str_replace_all(ALT_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = ifelse(ALT_AA == "", "del", ALT_AA))

  c = df %>%
    filter(EFFECT %in% c('conservative_inframe_insertion', 'disruptive_inframe_insertion')) %>%
    filter(grepl("dup", HGVS_P_modified)) %>%
    mutate(HGVS_P_modified = str_replace_all(HGVS_P_modified, setNames(replacement, patterns))) %>%
    separate(HGVS_P_modified, sep = "dup", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    mutate(ALT_AA = 'dup')

  d = df %>%
    filter(EFFECT == 'conservative_inframe_insertion') %>%
    filter(!grepl("dup", HGVS_P_modified)) %>%
    separate(HGVS_P_modified, sep = "ins", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    mutate(REF_AA = str_replace_all(REF_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = str_replace_all(ALT_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = paste("ins", ALT_AA, sep = ""))

  e1 = df %>%
    filter(EFFECT == 'disruptive_inframe_insertion') %>%
    filter(!grepl("dup", HGVS_P_modified)) %>%
    separate(HGVS_P_modified, sep = "del", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    filter(grepl("ins", ALT_AA)) %>%
    mutate(REF_AA = str_replace_all(REF_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = str_replace_all(ALT_AA, setNames(replacement, patterns)))

  e2 = df %>%
    filter(EFFECT == 'disruptive_inframe_insertion') %>%
    filter(!grepl("dup", HGVS_P_modified)) %>%
    filter(!grepl("del", HGVS_P_modified)) %>%
    separate(HGVS_P_modified, sep = "ins", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    mutate(REF_AA = str_replace_all(REF_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = str_replace_all(ALT_AA, setNames(replacement, patterns))) %>%
    mutate(ALT_AA = paste("ins", ALT_AA, sep = ""))

  f = df %>%
    filter(EFFECT == 'stop_gained') %>%
    separate(HGVS_P_modified, sep = "[0-9]+", into = c("REF_AA", "ALT_AA"), remove = F) %>%
    mutate(REF_AA = str_replace_all(REF_AA, setNames(replacement, patterns)))

  g = bind_rows(lapply(list(a, b, c, d, e1, e2, f), function(x) if(nrow(x) == 0) NULL else x)) %>%
    arrange(POS) %>%
    select(-HGVS_P_modified)

  return(g)
}



#' View HGVS_C Effects
#'
#' Plot the consequences of HGVS_C against the wildtype (NC_045512.2) gene sequence
#'
#' @param HGVS_C A string from snpSIFT output
#' @param gene The gene name
#' @param window An approximate upstream and downstream window size for plotting
#'
#' @examples
#' view_HGVS_C(HGVS_C = "c.326_328delCTT")
#'
#'
#' @export

view_HGVS_C = function(HGVS_C, gene = "S", window = 10) {
  require(tidyverse)
  require(Biostrings)
  require(msa)
  require(ggmsa)
  require(patchwork)

  genome_file = system.file("extdata", "NC_045512.2.gbk.ffn", package = "MappgeneSummary")
  seqs <- readDNAStringSet(genome_file, "fasta")

  # logic to do substitutions
  if (grepl("del", HGVS_C)) {
    n = stringr::str_extract_all(HGVS_C, "\\d+")
    start = as.numeric(n[[1]][1])
    stop = as.numeric(n[[1]][2])
    range = IRanges(start, stop)

    window_right = floor(start(range) / 3) * 3 + window * 3

    S = subseq(seqs[[gene]], 1,window_right)
    #S = seqs[[gene]]
    S_modified <- replaceAt(S, range, value=strrep("N", width(range)))

  } else if (grepl(">", HGVS_C)) {
    n = stringr::str_extract(HGVS_C, "\\d+")
    start = as.numeric(n[[1]][1])
    range = IRanges(start, start)

    window_right = floor(start(range) / 3) * 3 + window * 3

    S = subseq(seqs[[gene]], 1, window_right)
    #S = seqs[[gene]]
    S_modified <- replaceAt(S, range, value=str_sub(HGVS_C, -1))

  } else if (grepl("ins", HGVS_C)) {
    a = str_split(HGVS_C, pattern = "ins|dup")[[1]]
    insertion = a[2]
    n = stringr::str_extract_all(a[1], "\\d+")
    start = as.numeric(n[[1]][1])
    range = IRanges(start)
    y_range = IRanges(start, start+(nchar(insertion)))

    window_right = floor(start(range) / 3) * 3 + window * 3

    S = subseq(seqs[[gene]], 1, window_right)
    #S = seqs[[gene]]
    S <- replaceAt(S, range, value=paste0(as.character(subseq(S, start(range), end(range))), strrep("N", nchar(insertion))))
    S_modified <- replaceAt(S, y_range, value=paste0(as.character(subseq(S, start(range), end(range))), insertion))

    range = y_range # just to get the window size right

  } else if (grepl("dup", HGVS_C)) {
    a = str_split(HGVS_C, pattern = "ins|dup")[[1]]
    insertion = a[2]
    n = stringr::str_extract_all(a[1], "\\d+")
    start = as.numeric(n[[1]][2])
    range = IRanges(start)
    y_range = IRanges(start, start+(nchar(insertion)))

    window_right = floor(start(range) / 3) * 3 + window * 3

    S = subseq(seqs[[gene]], 1, window_right)
    #S = seqs[[gene]]
    S <- replaceAt(S, range, value=paste0(as.character(subseq(S, start(range), end(range))), strrep("N", nchar(insertion))))
    S_modified <- replaceAt(S, y_range, value=paste0(as.character(subseq(S, start(range), end(range))), insertion))

    range = y_range # just to get the window size right
  }

  aa_left = floor(start(range) / 3) - window
  if (aa_left < 1) {aa_left = 1}

  aa_right = ceiling(end(range) / 3) + window
  if (aa_right > length(seqs[[gene]])) {aa_right = length(seqs[[gene]])}

  gene_list <- list()
  gene_list[[paste(gene, "(NC_045512.2)")]] <- S
  gene_list[[paste(gene, "(modified)")]] <- S_modified

  dna = DNAStringSet(x = gene_list)
  aa = translate(dna, if.fuzzy.codon = "X")

  dna_aln = DNAMultipleAlignment(msa(dna, gapExtension = 0))
  aa_aln = AAMultipleAlignment(msa(aa, gapExtension = 0))


  a = ggmsa(dna, aa_left*3-2, aa_right*3, color = "Chemistry_NT", font = "DroidSansMono", char_width = 0.5, seq_name = TRUE) + geom_msaBar()
  b = ggmsa(aa_aln, aa_left, aa_right, color = "Chemistry_AA", font = "DroidSansMono", char_width = 0.5, seq_name = TRUE ) + geom_msaBar()

  a$plotlist[[2]] / b$plotlist[[2]] +
    plot_annotation(
      title = paste(gene, HGVS_C)
    )
}


