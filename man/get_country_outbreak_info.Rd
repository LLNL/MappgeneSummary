% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_outbreak_info.R
\name{get_country_outbreak_info}
\alias{get_country_outbreak_info}
\title{Country data from outbreak API}
\usage{
get_country_outbreak_info(countries, SNVs)
}
\arguments{
\item{countries}{A country or vector of countries}

\item{SNVs}{A SNV or vector of SNVs formatted for the outbreak API (eg S:Y1155H)}
}
\description{
The country code appears to be a three character string, like "USA" for the United States, or "CAN" for Canada.
}
\examples{
get_country_outbreak_info(countries = c("CA", "OR"), SNVs = "S:Y1155H")
}
