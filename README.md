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

## Customizing

There are several parameters that you can change in the .Rmd yaml header.

```
params:
  echo_type: FALSE
  mappgene_file_path: test
  summary_output_file_path: "." 
  overwrite_Rds_files: TRUE
  api_top_count: 10
```

`echo_type`: If `TRUE`, will show the code blocks in the html file.
`mappgene_file_path`: If set to 'test' it will use the test files included in the package. Otherwise, change this to the path where your mappgene output is located.
`summary_output_file_path`: The path location of the output files.
`overwrite_Rds_files`: Some code blocks take time to run. Set to `TRUE` to run all code blocks, or to `FALSE` to use saved Rds files.
`api_top_count`: Number of Outbreak API SNVs to generate plots for.
