#' Find the optimal min_score for searchSeq()
#'
#' \code{optimal_min_score()} finds the optimal min_score threshold for
#' TFBSTools::searchSeq by maximising the difference between the proportion
#' of motifs found in a subset of the peak data and a set of background control
#' sequences.
#'
#' @param peaks A GRanges peak object
#' @param pwm An object of class PWMatrix.
#' @param seed A single value specifying the state for random number generation.
#'
optimal_min_score <- function(peaks, pwm, genome_build, seed) {
  set.seed(seed)
  # subset peaks
  if (length(peaks) < 500) {
    subset_peaks <- peaks
  } else {
    subset_peaks <- sample(peaks, 500)
  }
  subset_ctrl <- adjacent_sequences(subset_peaks)

  loss_function <- function(param) {
    # calculate proportion for a given min_score
    calculate_score <- function(peaks,
                                pwm,
                                param,
                                genome_build) {
      hits <- get_motif_hits(
        peaks = peaks,
        pwm = pwm,
        min_score = param,
        genome_build = genome_build
      )
      sig <- TFBSTools::relScore(hits)

      hit_peak_names <- names(sig)[sapply(sig, length) > 0]
      hit_peaks <- peaks[names(peaks) %in% hit_peak_names,]

      return(length(hit_peaks) / length(peaks))
    }

    peak_score <- calculate_score(peaks = subset_peaks,
                                  pwm = pwm,
                                  param = param,
                                  genome_build)
    ctrl_score <- calculate_score(peaks = subset_ctrl,
                                  pwm = pwm,
                                  param = param,
                                  genome_build)

    return(ctrl_score - peak_score)  # Invert to find maximum
  }

  result <- stats::optimise(loss_function,
                            lower = 0,
                            upper = 1,
                            maximum = FALSE)

  return(result$minimum)
}
