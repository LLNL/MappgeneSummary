# library(devtools)
# library(roxygen2)

devtools::document()

devtools::install("../MappgeneSummary")
packageVersion("MappgeneSummary")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("Biostrings")
#  BiocManager::install("msa")
#  BiocManager::install("ggmsa")

#from github
#install_github("jeffkimbrel/MappgeneSummary")
install_github("jeffkimbrel/jakR")

MappgeneSummary::copy_summary_template("~/Desktop/", name = "test.Rmd")
file.edit("~/Desktop/test.Rmd")
#MappgeneSummary::get_state_outbreak_info(states = c("ca", "OR"), SNVs = c("S:y1155H", "S:P681H"))
#MappgeneSummary::get_country_outbreak_info(countries = c("usa", "can"), SNVs = c("S:y1155H", "S:P681H"))
#MappgeneSummary::get_global_outbreak_info(SNVs = c("S:y1155H", "S:P681H"))


# New outbreak.info API and authentication
#library(outbreakinfo)
#authenticateUser()

# vignette (https://r-pkgs.org/vignettes.html)
usethis::use_vignette("my-vignette")
