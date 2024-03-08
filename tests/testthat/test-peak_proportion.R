test_that("list outputted by peak_proportion()", {
  peak_file <- system.file("extdata",
                           "rep1_peaks.narrowPeak",
                           package = "MotifStats")
  data("creb_motif", package = "MotifStats")

  res <- peak_proportion(
    peak_file = peak_file,
    pwm = creb_motif,
    plot = FALSE,
    min_score = 0.8
  )

  expect_true(is.list(res))
})

ctcf_peaks <- "read1_no_ctrl_01_peaks.narrowPeak"
motif_pwm <- read_motif_file(motif_file = "MA1930.2.jaspar",
                             file_format = "jaspar")

peak_proportion(peak_file = ctcf_peaks,
                pwm = motif_pwm,
                min_score = 0.6)



