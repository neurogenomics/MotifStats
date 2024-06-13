test_that("The output of read_peak_file is a GRanges object", {
  peak_file <- system.file("extdata",
                           "ctcf_seacr.stringent.bed",
                           package = "MotifStats")
  res <- read_peak_file_seacr(file_path = peak_file)

  expect_s4_class(res, class = "GRanges")
})
