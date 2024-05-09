test_that("the density_plot function outputs a ggplot", {
  temp_dir <- withr::local_tempdir()

  data("ctcf_peaks", package = "MotifStats")
  data("ctcf_motif", package = "MotifStats")

  summit_to_motif_out <- summit_to_motif(
    peak_input = ctcf_peaks,
    motif = ctcf_motif,
    fp_rate = 0.05,
    out_dir = temp_dir,
    genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )
  res <- density_plot(
    distance_vec = summit_to_motif_out$distance_to_summit)

  expect_true(is.ggplot(res))
})
