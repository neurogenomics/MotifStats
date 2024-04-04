test_that("The output of read_peak_file is a GRanges object", {
  peak_file <- system.file("extdata",
                           "rep1_peaks.narrowPeak",
                           package = "MotifStats")
  res <- read_peak_file(file_path = peak_file)

  expect_s4_class(res, class = "GRanges")
})
