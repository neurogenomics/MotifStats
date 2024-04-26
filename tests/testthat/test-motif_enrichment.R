test_that("list outputted by peak_proportion function", {
  data("creb_peaks", package = "MotifStats")
  data("creb_motif", package = "MotifStats")

  temp_dir <- withr::local_tempdir()
  withr::defer(temp_dir)

  res <- motif_enrichment(
    peak_input = creb_peaks,
    motif = creb_motif,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    out_dir = temp_dir
    )
  expect_true(is.list(res))
})
