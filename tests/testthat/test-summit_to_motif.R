test_that("list outputted by summit_to_motif function", {
  temp_dir <- withr::local_tempdir()

  data("ctcf_peaks", package = "MotifStats")
  data("ctcf_motif", package = "MotifStats")

  res <- summit_to_motif(
    peak_input = ctcf_peaks,
    motif = ctcf_motif,
    fp_rate = 5e-02,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    out_dir = temp_dir
    )
  expect_true(is.list(res))
})
