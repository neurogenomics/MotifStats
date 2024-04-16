creb_motif <- read_motif_file(
  "inst/extdata/MA0018.5.jaspar",
  motif_id = "MA0018.5",
  file_format = "jaspar"
  )
usethis::use_data(creb_motif, overwrite = TRUE)
