library("devtools")
library("roxygen2")

document()

install("../MappgeneSummary")
packageVersion("MappgeneSummary")

#from github
#install_github("jeffkimbrel/MappgeneSummary")
install_github("jeffkimbrel/jakR")

MappgeneSummary::copy_summary_template("~/Desktop/", name = "test.Rmd")
