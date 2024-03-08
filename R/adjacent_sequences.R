#' Generate sequences immediately adjacent to peaks
#'
#' \code{adjacent_sequences()} generates sequences adjacent to each peak in
#' "peaks". These adjacent sequences will have the same width as the peaks.
#' These are used as a reference when calculating the proportion of peaks with
#' a motif. Where a peak is at the boundaries of a chromosome the adjacent
#' sequence width will not go beyond the boundaries of the chromosome.
#'
#' @param peaks A GRanges peak object
#'
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlengths
#'
#' @keywords internal
adjacent_sequences <- function(peaks) {
  # Retrieve the width of each peak
  peak_widths <- GenomicRanges::width(peaks)

  # Retrieve chromosome lengths for boundary checks
  chrom_sizes <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)

  # Randomly choose the side that each adjacent sequence will be placed
  is_upstream <- sample(c(TRUE, FALSE), size = length(peaks), replace = TRUE)

  for (i in seq_along(peaks)) {
    chrom <- as.character(seqnames(peaks)[i])
    chrom_length <- chrom_sizes[chrom]

    if (is_upstream[i]) {
      # Control sequence will be upstream
      new_end <- start(peaks)[i] - 1
      new_start <- max(1, new_end - peak_widths[i] + 1)
    } else {
      # Control sequence will be downstream
      new_start <- end(peaks)[i] + 1
      new_end <- min(chrom_length, new_start + peak_widths[i] - 1)
    }

    # Update the peak coordinates
    start(peaks)[i] <- new_start
    end(peaks)[i] <- new_end
  }

  return(peaks)
}
