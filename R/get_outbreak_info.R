#' Data from outbreak API
#'
#' @param SNVs A SNV or vector of SNVs formatted for the outbreak API (eg S:Y1155H)
#'
#' @export
#'
#' @examples
#' get_outbreak_info(SNVs = "S:Y1155H")
#'

get_outbreak_info = function(SNVs) {
  df = data.frame()

  for (SNV in SNVs) {
    California = getPrevalence(mutations = SNV, location = "California", logInfo = FALSE)
    USA = getPrevalence(mutations = SNV, location = "United States", logInfo = FALSE)
    Global = getPrevalence(mutations = SNV, logInfo = FALSE)

    df = rbind(df, California, USA, Global)
  }

  return(df)
}

#' State data from outbreak API
#'
#' @param states A state or vector of states
#' @param SNVs A SNV or vector of SNVs formatted for the outbreak API (eg S:Y1155H)
#'
#' @export
#'
#' @examples
#' get_state_outbreak_info(states = c("CA", "OR"), SNVs = "S:Y1155H")

get_state_outbreak_info = function(states, SNVs) {

  require(tidyverse)
  require(jsonlite)

  df = data.frame()

  for (state in states) {
    state = toupper(state)

    for (SNV in SNVs) {
      SNV = toupper(SNV)

      state_data = tryCatch(data.frame(fromJSON(paste0("https://api.outbreak.info/genomics/prevalence-by-location?location_id=USA_US-", state, "&mutations=", SNV))),
                            error = function(e)
                            data.frame("success" = FALSE))

      state_data = data.frame(state_data) %>%
        mutate("LOCATION" = state)

      names(state_data) = gsub(x = names(state_data), pattern = paste0("results.", gsub(x = SNV, pattern = ":", replacement = "."), "."), replacement = "")

      if (!'date' %in% colnames(state_data)) {
        message(paste(SNV, 'not found in', state))
      } else {
        state_data = state_data %>%
          mutate(date = as.Date(date), SNV = SNV)
      }

      df = bind_rows(df, state_data)
    }
  }

  return(as_tibble(df))
}

#' Country data from outbreak API
#'
#' The country code appears to be a three character string, like "USA" for the United States, or "CAN" for Canada.
#'
#' @param countries A country or vector of countries
#' @param SNVs A SNV or vector of SNVs formatted for the outbreak API (eg S:Y1155H)
#'
#' @export
#'
#' @examples
#' get_country_outbreak_info(countries = c("CA", "OR"), SNVs = "S:Y1155H")

get_country_outbreak_info = function(countries, SNVs) {

  require(tidyverse)
  require(jsonlite)

  df = data.frame()

  for (country in countries) {
    country = toupper(country)

    for (SNV in SNVs) {
      SNV = toupper(SNV)

      country_data = tryCatch(data.frame(fromJSON(paste0("https://api.outbreak.info/genomics/prevalence-by-location?location_id=", country, "&mutations=", SNV))),
                            error = function(e)
                              data.frame("success" = FALSE))

      country_data = data.frame(country_data) %>%
        mutate("LOCATION" = country)

      names(country_data) = gsub(x = names(country_data), pattern = paste0("results.", gsub(x = SNV, pattern = ":", replacement = "."), "."), replacement = "")

      if (!'date' %in% colnames(country_data)) {
        message(paste(SNV, 'not found in', country))
      } else {
        country_data = country_data %>%
          mutate(date = as.Date(date), SNV = SNV)
      }

      df = bind_rows(df, country_data)
    }
  }

  return(as_tibble(df))
}

#' Global data from outbreak API
#'
#' @param SNVs A SNV or vector of SNVs formatted for the outbreak API (eg S:Y1155H)
#'
#' @export
#'
#' @examples
#' get_global_outbreak_info(SNVs = "S:Y1155H")

get_global_outbreak_info = function(SNVs) {

  require(tidyverse)
  require(jsonlite)

  df = data.frame()

  for (SNV in SNVs) {
    SNV = toupper(SNV)

    global_data = tryCatch(data.frame(fromJSON(paste0("https://api.outbreak.info/genomics/prevalence-by-location?mutations=", SNV))),
                            error = function(e)
                              data.frame("success" = FALSE))

    global_data = data.frame(global_data) %>%
      mutate("LOCATION" = "GLOBAL")

    names(global_data) = gsub(x = names(global_data), pattern = paste0("results.", gsub(x = SNV, pattern = ":", replacement = "."), "."), replacement = "")

    if (!'date' %in% colnames(global_data)) {
      message(paste(SNV, 'not found in GLOBAL'))
    } else {
      global_data = global_data %>%
        mutate(date = as.Date(date), SNV = SNV)
    }

    df = bind_rows(df, global_data)

  }

  return(as_tibble(df))
}
