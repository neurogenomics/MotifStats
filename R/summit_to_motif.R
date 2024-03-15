#' Calculate the distance between peak summits and motifs
#'
#' \code{summit_to_motif()} calculates the distance between each motif and its
#' nearest summit. A peak summit represents the exact nucleotide where signal
#' enrichment is highest.
#'
#' @import rtracklayer
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import TFBSTools
#'
#' @inheritParams peak_proportion
#'
#' @returns GRanges peak object with metadata columns relating to motif
#' positions.
#'
#' @examples
#' \dontrun{
#' pwm <-
#'   read_motif_file(motif_file = "./MA0018.3.meme", file_format = "meme")
#'
#' summit_to_motif(peak_file = "read1_01_peaks.narrowPeak",
#'                 pwm = pwm,
#'                 min_score = NULL,
#'                 optimal_min_score = TRUE,
#'                 seed = 123
#'                 )
#'                 }
#' @export
summit_to_motif <- function(peak_input,
                            pwm,
                            min_score = 0.8,
                            genome_build,
                            optimal_min_score = FALSE,
                            seed = NULL) {
  if(!is.null(seed)) set.seed(seed)

  if (is.character(peak_input)) {
    normalizePath(peak_input, mustWork = "TRUE") # error if path is invalid
    peaks <- read_peak_file(peak_input)
  } else if (inherits(peak_input, "GRanges")) {
    peaks <- peak_input
  } else {
    stopper(
      "peak_input must be a path to a peak file or a GRanges object",
      "with a 0-based summit column. This can be generated using",
      "read_peak_file(file_path)."
    )
  }

  if (optimal_min_score) {
    min_score <- optimal_min_score(
      peaks = peaks,
      pwm = pwm,
      genome_build = genome_build
    )
    messager("Optimal min_score found:", min_score)

  } else if (!is.null(min_score)) {
    min_score <- min_score

  } else {
    stopper("Either optimal_min_score is set to TRUE or min_score must be",
            "supplied.")
  }

  hits <- get_motif_hits(
    peaks = peaks,
    pwm = pwm,
    min_score = min_score,
    genome_build = genome_build
  )
  # Write the motif hits to GFF3
  motif_coord_df <- TFBSTools::writeGFF3(hits)

  # Find peaks with multiple motifs
  index_to_repeat <-
    base::match(motif_coord_df$seqname, names(peaks))

  # Expand the peak object to repeat the peaks that have >1 motif
  expanded_peaks <- peaks[index_to_repeat]

  # Add the motif start and end data for each peak
  mcols(expanded_peaks)$motif_start <-
    GenomicRanges::start(expanded_peaks) + motif_coord_df$start
  mcols(expanded_peaks)$motif_end <-
    GenomicRanges::start(expanded_peaks) + motif_coord_df$end

  # Calculate motif centres and distance to summit
  motif_center <-
    (mcols(expanded_peaks)$motif_start + mcols(expanded_peaks)$motif_end) / 2
  distance_to_summit <-
    abs(mcols(expanded_peaks)$summit - motif_center)

  # Adding distance to metadata
  mcols(expanded_peaks)$distance_to_summit <- distance_to_summit

  return(
    list(
      peak_set = expanded_peaks,
      distance_to_summit = distance_to_summit,
      seed = seed,
      min_score = min_score
    )
  )
}

# CTCF motif
# pwm <- read_motif_file("MA1930.2.jaspar", file_format = "jaspar")
# one <- summit_to_motif("./both_01_peaks.narrowPeak",
#                         pwm = pwm)
# two <- summit_to_motif("./read1_01_peaks.narrowPeak",
#                         pwm = pwm)
# par(mfrow = c(1,2))
# hist(one$distance_to_summit, breaks = 40, xlim = c(0, 1000))
# hist(two$distance_to_summit, breaks = 40, xlim = c(0, 1000))
#
# median(one$distance_to_summit)
# median(two$distance_to_summit)
