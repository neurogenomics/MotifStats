#' Generate sequences immediately adjacent to peaks
#'
#' \code{adjacent_sequences()} generates sequences adjacent to each peak in
#' "peaks".
#'
#' Adjacent sequences will have the same width as their respective
#' peaks. Where a peak is at the boundary of a chromosome, the adjacent
#' sequence width will not extend beyond the boundaries of the chromosome.
#'
#' @param peaks A GRanges peak object
#'
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomeInfoDb seqlengths
#'
#' @keywords internal

adjacent_sequences <- function(peaks,
                               genome_build) {
  peak_widths <- GenomicRanges::width(peaks)
  chrom_sizes <- GenomeInfoDb::seqlengths(genome_build)
  is_upstream <-
    sample(c(TRUE, FALSE), size = length(peaks), replace = TRUE)

  # Define a function to process each peak
  process_peak <- function(i) {
    chrom <- as.character(seqnames(peaks)[i])
    chrom_length <- chrom_sizes[chrom]
    is_upstream <- is_upstream[i]
    peak_width <- peak_widths[i]

    if (is_upstream) {
      # adjacent sequence will be upstream
      new_end <- start(peaks)[i] - 1
      new_start <- new_end - peak_width + 1
      if (new_start < 1) {
        # if the start is out of bounds (try downstream)
        new_start <- end(peaks)[i] + 1
        new_end <- new_start + peak_width - 1
        if (new_end > chrom_length) {
          # if the end is out of bounds
          new_start <- end(peaks)[i] + 1
          new_end <- min(chrom_length, new_start + peak_width - 1)
        }
      }
    } else {
      # adjacent sequence will be downstream
      new_start <- end(peaks)[i] + 1
      new_end <- new_start + peak_width - 1
      if (new_end > chrom_length) {
        # if the end is out of bounds (try upstream)
        new_end <- start(peaks)[i] - 1
        new_start <- new_end - peak_width + 1
        if (new_start < 1) {
          # if the start is out of bounds
          new_end <- start(peaks)[i] - 1
          new_start <- max(new_end - peak_width + 1, 1)
        }
      }
    }

    return(c(new_start, new_end))
  }
  # Apply the function to each peak
  peak_background <- t(sapply(seq_along(peaks), process_peak))

  # Update peak positions
  start(peaks) <- peak_background[, 1]
  end(peaks) <- peak_background[, 2]

  return(peaks)
}
