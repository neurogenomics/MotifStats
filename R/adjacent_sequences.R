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

# adjacent_sequences <- function(peaks) {
#   # Retrieve the width of each peak using width() accessor
#   peak_widths <- GenomicRanges::width(peaks)
#
#   # Retrieve chromosome lengths for boundary checks
#   chrom_sizes <-
#     GenomeInfoDb::seqlengths( # add genome_build option
#       BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#       )
#
#   # Randomly choose the side that each adjacent sequence will be placed
#   is_upstream <- sample(c(TRUE, FALSE), size = length(peaks), replace = TRUE)
#
#   for (i in seq_along(peaks)) {
#     chrom <- as.character(seqnames(peaks)[i])
#     chrom_length <- chrom_sizes[chrom]
#
#     if (is_upstream[i]) {
#       # Control sequence will be upstream
#       new_end <- start(peaks)[i] - 1
#       new_start <- max(1, new_end - peak_widths[i] + 1)
#     } else {
#       # Control sequence will be downstream
#       new_start <- end(peaks)[i] + 1
#       new_end <- min(chrom_length, new_start + peak_widths[i] - 1)
#     }
#
#     # Update the peak coordinates
#     start(peaks)[i] <- new_start
#     end(peaks)[i] <- new_end
#   }
#
#   return(peaks)
# }

adjacent_sequences <- function(peaks) {
  peak_widths <- GenomicRanges::width(peaks)
  chrom_sizes <- GenomeInfoDb::seqlengths(
    BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )
  is_upstream <- sample(c(TRUE, FALSE), size = length(peaks), replace = TRUE)

  # Define a function to process each peak
  process_peak <- function(i) {
    chrom <- as.character(seqnames(peaks)[i])
    chrom_length <- chrom_sizes[chrom]
    is_upstream <- is_upstream[i]
    peak_width <- peak_widths[i]

    if (is_upstream) {
      new_end <- start(peaks)[i] - 1
      new_start <- new_end - peak_width + 1
      if (new_start < 1) { # If the position of the seq start is < 1 (out of bounds)
        new_start <- end(peaks)[i] + 1
        new_end <- new_start + peak_width - 1
      }
    } else {
      new_start <- end(peaks)[i] + 1
      new_end <- new_start + peak_width - 1
      if (new_end > chrom_length) { # if the position of the seq end is > the chrom length (out of bounds)
        new_end <- start(peaks)[i] - 1
        new_start <- new_end - peak_width + 1
      }
    }

    c(new_start, new_end)
  }

  # Apply the function to each peak
  peak_background <- t(sapply(seq_along(peaks), process_peak))

  # Update peak positions
  start(peaks) <- peak_background[, 1]
  end(peaks) <- peak_background[, 2]

  return(peaks)
}
