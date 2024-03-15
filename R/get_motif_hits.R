#' Scan a set of peaks for a specified motif
#'
#' \code{get_motif_hits()} uses the \code{searchSeq()} function from
#' \code{TFBSTools} to scan each peak for the presence of a motif.
#'
#' @import BSgenome
#' @importFrom TFBSTools searchSeq
#'
#' @inheritParams peak_proportion
#' @param peaks A GRanges peak object.
#'
#' @keywords internal
get_motif_hits <- function(peaks,
                           pwm,
                           min_score,
                           genome_build = genome_build){
  peak_sequences <-
    BSgenome::getSeq(genome_build,
                     peaks)

  hits <- TFBSTools::searchSeq( # Consider using scan_sequence instead of searchSeq.
    x = pwm,
    subject = peak_sequences,
    strand = "*",
    min.score = min_score
  )

  return(hits)
}
