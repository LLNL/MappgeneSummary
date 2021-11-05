library("devtools")
library("roxygen2")

document()

install("../MappgeneSummary")
packageVersion("MappgeneSummary")

#from github
#install_github("jeffkimbrel/MappgeneSummary")
#install_github("jeffkimbrel/jakR")

#MappgeneSummary::copy_summary_template("~/Desktop/", name = "test.Rmd")
#MappgeneSummary::get_state_outbreak_info(states = c("ca", "OR"), SNVs = c("S:y1155H", "S:P681H"))
#MappgeneSummary::get_country_outbreak_info(countries = c("usa", "can"), SNVs = c("S:y1155H", "S:P681H"))
#MappgeneSummary::get_global_outbreak_info(SNVs = c("S:y1155H", "S:P681H"))
