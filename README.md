# MappgeneSummary

# Installation from Github

```
install_github("jeffkimbrel/MappgeneSummary")
```

Additionally, there are a few functions and color palettes from the `jakR` package. Install from github...

```
install_github("jeffkimbrel/jakR")
```

# Examples

```
get_state_outbreak_info(states = c("ca", "OR"), SNVs = c("S:y1155H", "S:P681H"))
get_country_outbreak_info(countries = c("usa", "can"), SNVs = c("S:y1155H", "S:P681H"))
get_global_outbreak_info(SNVs = c("S:y1155H", "S:P681H"))
```

# Summary.Rmd File

To make a copy of the `.Rmd` file and the required `styles.css` file, type the following giving a path and file name. By default, the html and file outputs will be written to the same `path`. 

```
MappgeneSummary::copy_summary_template(path = "~/Desktop/", name = "test.Rmd")
```
