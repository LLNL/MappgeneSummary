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
