ctcf_motif <- read_motif_file(
  "inst/extdata/MA1930.2.jaspar",
  motif_id = "MA1930.2.jaspar",
  file_format = "jaspar"
)
usethis::use_data(ctcf_motif, overwrite = TRUE)
