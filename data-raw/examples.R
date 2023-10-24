## code to prepare `examples` dataset goes here

example_df_A_lofreq = readr::read_tsv("inst/test_files/sampleA.ivar.lofreq.snpSIFT.txt")
example_df_B_lofreq = readr::read_tsv("inst/test_files/sampleB.ivar.lofreq.snpSIFT.txt")
example_df_C_lofreq = readr::read_tsv("inst/test_files/sampleC.ivar.lofreq.snpSIFT.txt")
example_df_A_ivar = readr::read_tsv("inst/test_files/sampleA.ivar.snpSIFT.txt")
example_df_B_ivar = readr::read_tsv("inst/test_files/sampleB.ivar.snpSIFT.txt")
example_df_C_ivar = readr::read_tsv("inst/test_files/sampleC.ivar.snpSIFT.txt")

usethis::use_data(example_df_A_lofreq,
                  example_df_B_lofreq,
                  example_df_C_lofreq,
                  example_df_A_ivar,
                  example_df_B_ivar,
                  example_df_C_ivar,
                  overwrite = TRUE)
