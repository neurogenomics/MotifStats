test_that("list outputted by motif_enrichment function", {
  data("ctcf_peaks", package = "MotifStats")
  data("ctcf_motif", package = "MotifStats")

  temp_dir <- withr::local_tempdir()
  withr::defer(temp_dir)

  res <- motif_enrichment(
    peak_input = ctcf_peaks,
    motif = ctcf_motif,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    out_dir = temp_dir
    )
  expect_true(is.list(res))
})
