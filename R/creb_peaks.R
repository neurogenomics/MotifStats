#' CREB Peaks
#'
#' This dataset contains a set of CREB1 TIP-seq peaks (narrowPeak) produced by
#' MACS3 (first-mate, 5' shift). We have subset the peaks to reduce the data
#' size (only chr19). The commands used to subset the data were:
#'
#' creb_peaks <- read.table("path/to/creb_peaks")
#' creb_peaks <- creb_peaks[creb_peaks$V1 == "chr19",]
#' write.table(creb_peaks,
#'             "creb_subset.narrowPeak",
#'             row.names = FALSE,
#'             col.names = FALSE,
#'             quote = FALSE,
#'             sep = "\t")
#'
#' @format A GRanges peak object outputted by the \code{read_peak_file}
#' function.
"creb_peaks"
