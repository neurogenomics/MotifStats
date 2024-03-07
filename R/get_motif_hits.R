#' Scan a set of peaks for a specified motif
#'
#' \code{get_motif_hits()} uses the \code{searchSeq()} function from
#' \code{TFBSTools} to scan each peak for the presence of a motif.
#'
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom TFBSTools searchSeq
#'
#' @inheritParams peak_proportion
#' @param peaks A GRanges peak object.
#'
#' @keywords internal
get_motif_hits <- function(peaks,
                           pwm,
                           min_score){
  peak_sequences <-
    BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                     peaks)
  hits <- TFBSTools::searchSeq(
    x = pwm,
    subject = peak_sequences,
    strand = "*",
    min.score = min_score
  )

  return(hits)
}
