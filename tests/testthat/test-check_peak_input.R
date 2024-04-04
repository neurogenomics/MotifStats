test_that("fail when peaks are incorrectly formatted", {
  not_a_peak <- data.frame(a = c(1:10),
                           b = runif(10),
                           c = rnorm(10))
  expect_error(check_peak_input(
    not_a_peak, # expects narrowPeak or GRanges
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ))
})

test_that("list outputted by check_peak_input function", {
  data("creb_peaks", package = "MotifStats")
  res <- check_peak_input(
    creb_peaks,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )

  expect_true(is.list(res))
})
