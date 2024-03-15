test_that("list outputted by peak_proportion()", {
  peak_file <- system.file("extdata",
                           "rep1_peaks.narrowPeak",
                           package = "MotifStats")
  data("creb_motif", package = "MotifStats")

  res <- peak_proportion(
    peak_input = peak_file,
    pwm = creb_motif,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    min_score = 0.8,
    optimal_min_score = FALSE,
    seed = 123
  )

  expect_true(is.list(res))
})
