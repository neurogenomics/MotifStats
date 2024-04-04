test_that("the output of read_motif_file is a universalmotif object", {
  motif_file <- system.file("extdata",
                           "MA0018.5.jaspar",
                           package = "MotifStats")
  res <- read_motif_file(motif_file = motif_file,
                         motif_id = "CREB1",
                         file_format = "jaspar")

  expect_s4_class(res, class = "universalmotif")
})

test_that("fail when unsupported file format is given", {
  motif_file <- system.file("extdata",
                            "MA0018.5.jaspar",
                            package = "MotifStats")
  expect_error(read_motif_file(
    motif_file = motif_file,
    motif_id = "CREB1",
    file_format = "invalid_format"
  ))
})
