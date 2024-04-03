#' Check that the peak input is valid
#'
#' \code{check_peak_input} confirms that the peak input is valid.
#'
#' @importFrom BSgenome getSeq
#'
#' @inheritParams peak_proportion
#'
#' @returns A list containing a GRanges peak object and a DNAStringSet sequence
#' object.
#'
#' @keywords internal

check_peak_input <- function(peak_input,
                             genome_build){

  if (is.character(peak_input)) { # If peak_input is a path to a peak file
    normalizePath(peak_input, mustWork = "TRUE")
    peaks <- read_peak_file(peak_input)
  } else if (inherits(peak_input, "GRanges")) {
    peaks <- peak_input
  } else {
    stopper(
      "peak_input must be a path to a peak file or a GRanges object",
      "This can be generated using read_peak_file(file_path)."
    )
  }
  peak_sequences <- BSgenome::getSeq(genome_build, peaks)

  return(list(peaks, peak_sequences))
}
