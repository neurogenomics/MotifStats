#' Generate sequences immediately adjacent to peaks
#'
#' \code{adjacent_sequences()} generates sequences adjacent to each peak in
#' "peaks". These sequences will be used as a reference when calculating the
#' proportion of peaks with a motif.
#'
#' @param peaks A GRanges peak object
#'
#' @keywords internal
adjacent_sequences <- function(peaks){
  # Calculate the width of each peak
  peak_widths <- GenomicRanges::width(peaks)

  # Randomly choose direction for each peak: TRUE for upstream, FALSE for downstream
  is_upstream <- sample(c(TRUE, FALSE), size = length(peaks), replace = TRUE)

  # Update starts and ends based on direction
  for (i in seq_along(peaks)) {
    if (is_upstream[i]) {
      # control sequence will be upstream
      new_end <- start(peaks)[i] - 1
      new_start <- new_end - peak_widths[i] + 1
    } else {
      # controls sequence will be downstream
      new_start <- end(peaks)[i] + 1
      new_end <- new_start + peak_widths[i] - 1
    }

    # Update the peak coordinates
    start(peaks)[i] <- new_start
    end(peaks)[i] <- new_end
  }

  return(peaks)
}
